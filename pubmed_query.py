"""Batch PubMed queries by author name or PMID using NCBI E-utilities (Biopython)."""

import os
import time
from urllib.error import HTTPError

from Bio import Entrez

Entrez.email = os.environ.get("NCBI_EMAIL", "you@example.com")  # required by NCBI
Entrez.api_key = os.environ.get("NCBI_API_KEY")  # optional, raises limit 3 -> 10 req/s

BATCH_SIZE = 200
DELAY = 0.11 if Entrez.api_key else 0.34  # stay under NCBI rate limits


def _call(func, retries=3, **kwargs):
    """Call an Entrez function, parse the result, and retry on transient errors."""
    for attempt in range(retries):
        try:
            with func(**kwargs) as handle:
                result = Entrez.read(handle)
            time.sleep(DELAY)
            return result
        except (HTTPError, RuntimeError) as e:
            if attempt == retries - 1:
                raise
            print(f"Entrez error ({e}); retrying...")
            time.sleep(2 ** attempt)


def _parse_article(article):
    """Flatten one PubmedArticle record into a dict."""
    citation = article["MedlineCitation"]
    art = citation["Article"]

    authors = []
    for a in art.get("AuthorList", []):
        if "CollectiveName" in a:
            authors.append(str(a["CollectiveName"]))
        else:
            authors.append(f"{a.get('LastName', '')} {a.get('Initials', '')}".strip())

    abstract_parts = []
    for part in art.get("Abstract", {}).get("AbstractText", []):
        label = getattr(part, "attributes", {}).get("Label")
        abstract_parts.append(f"{label}: {part}" if label else str(part))

    pub_date = art["Journal"]["JournalIssue"]["PubDate"]
    year = pub_date.get("Year") or pub_date.get("MedlineDate", "")[:4]

    doi = next(
        (str(x) for x in article["PubmedData"].get("ArticleIdList", [])
         if x.attributes.get("IdType") == "doi"),
        None,
    )

    return {
        "pmid": str(citation["PMID"]),
        "title": str(art.get("ArticleTitle", "")),
        "abstract": " ".join(abstract_parts),
        "authors": authors,
        "journal": str(art["Journal"].get("Title", "")),
        "year": str(year),
        "doi": doi,
        "mesh_terms": [str(m["DescriptorName"]) for m in citation.get("MeshHeadingList", [])],
        "pub_types": [str(p) for p in art.get("PublicationTypeList", [])],
    }


def fetch_pmids(pmids):
    """Fetch full records for a list of PMIDs, in batches. Returns a list of dicts."""
    pmids = list(dict.fromkeys(str(p).strip() for p in pmids if str(p).strip()))
    records = []
    for i in range(0, len(pmids), BATCH_SIZE):
        chunk = pmids[i:i + BATCH_SIZE]
        result = _call(Entrez.efetch, db="pubmed", id=",".join(chunk),
                       rettype="xml", retmode="xml")
        records.extend(_parse_article(a) for a in result["PubmedArticle"])
        print(f"Fetched {min(i + BATCH_SIZE, len(pmids))}/{len(pmids)} PMIDs")
    return records


def search_names(names, extra_query=None, max_per_name=500):
    """Search PubMed by author name; returns {name: [pmids]}.

    Names work best as "Lastname Initials" (e.g. "Smith JA"). `extra_query`
    is ANDed onto every search, e.g. '2020:2026[dp] AND cancer[mh]'.
    """
    hits = {}
    for name in names:
        term = f'"{name}"[au]'
        if extra_query:
            term += f" AND ({extra_query})"
        result = _call(Entrez.esearch, db="pubmed", term=term, retmax=max_per_name)
        hits[name] = list(result["IdList"])
        total = int(result["Count"])
        note = f" (capped at {max_per_name})" if total > max_per_name else ""
        print(f"{name}: {total} hits{note}")
    return hits


def batch_query(items, mode="auto", extra_query=None, max_per_name=500):
    """Query PubMed with a list of PMIDs or author names.

    mode: "pmid", "name", or "auto" (PMIDs if every item is numeric).
    Returns a list of article dicts. In name mode each dict also has
    "matched_names" listing which input names it came from.
    """
    items = [str(x).strip() for x in items if str(x).strip()]
    if mode == "auto":
        mode = "pmid" if all(x.isdigit() for x in items) else "name"

    if mode == "pmid":
        return fetch_pmids(items)
    if mode != "name":
        raise ValueError(f"mode must be 'pmid', 'name', or 'auto', got {mode!r}")

    hits = search_names(items, extra_query, max_per_name)
    pmid_to_names = {}
    for name, pmids in hits.items():
        for pmid in pmids:
            pmid_to_names.setdefault(pmid, []).append(name)

    records = fetch_pmids(pmid_to_names)
    for r in records:
        r["matched_names"] = pmid_to_names.get(r["pmid"], [])
    return records


if __name__ == "__main__":
    import pandas as pd

    df = pd.DataFrame(batch_query(["31452104", "33301246"]))
    print(df[["pmid", "year", "journal", "title"]])

    df = pd.DataFrame(batch_query(["Topol EJ"], extra_query="2024[dp]", max_per_name=5))
    print(df[["pmid", "matched_names", "title"]])

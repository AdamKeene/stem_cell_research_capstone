# Stem Cell Research Capstone

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Set your NCBI credentials (an API key is optional, but it raises the rate limit from 3 to 10 requests per second):

```bash
export NCBI_EMAIL="you@example.com"
export NCBI_API_KEY="your_key"
```

## PubMed queries

`pubmed_query.py` batch-queries PubMed by PMID or author name:

```python
from pubmed_query import batch_query
import pandas as pd

df = pd.DataFrame(batch_query(["31452104", "33301246"]))               # PMIDs
df = pd.DataFrame(batch_query(["Smith JA"], extra_query="2020:2026[dp]"))  # author names
```

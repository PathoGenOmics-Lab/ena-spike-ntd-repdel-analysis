import sys
import csv
from typing import Mapping, List, Iterator
import json
import requests
import urllib.parse

import pandas as pd


SEARCH_FIELDS_URL = "https://www.ebi.ac.uk/ena/portal/api/returnFields?dataPortal=ena&result=read_run&format=json"
POST_URL = "https://www.ebi.ac.uk/ena/portal/api/search"

HEADERS = {
    "content-type" : "application/x-www-form-urlencoded",
    "connection": "keep-alive"
}
REQUEST_DATA_TEMPLATE = {
    "limit": snakemake.params.limit,
    "result": "read_run",
    "format": "tsv"
}


def get_ENA_text_fields():
    r = requests.get(SEARCH_FIELDS_URL)
    if r.status_code != 200:
        sys.exit(f"<warning> Could not get search fields from ENA (code: {r.status_code})")
    fields = []
    for field in json.loads(r.text):
        if (name := field.get("columnId", None)): # and (field.get("type", None) == "text"
            fields.append(name)
    return fields


def build_query(min_date: str, max_date: str, tax_id: str, host: str) -> str:
    items = [
        f"tax_tree({tax_id})" if tax_id else None,
        f"collection_date>={min_date}" if min_date else None,
        f"collection_date<={max_date}" if max_date else None,
        f'host_scientific_name="{host}"' if host else None
    ]
    return " AND ".join(item for item in items if item)


def write_iter(content: Iterator, path: str):
    with open(path, "wb") as fw:
        for chunk in content:
            fw.write(chunk)


if __name__ == "__main__":

    # Init search
    print("Searching ENA fields", flush=True)
    fields = get_ENA_text_fields()

    # Perform search
    print("Performing ENA search", flush=True)
    data = REQUEST_DATA_TEMPLATE.copy()
    query_txt = build_query(
        snakemake.params.start_date,
        snakemake.params.end_date,
        snakemake.params.taxonomy,
        snakemake.params.host_scientific_name
    )
    data["query"] = urllib.parse.quote(query_txt, safe="")
    data["fields"] = ",".join(fields)
    with requests.post(POST_URL, data=data, headers=HEADERS, stream=True) as r:
        if r.status_code == 200:
            print("Writing ENA search results", flush=True)
            write_iter(r.iter_content(chunk_size=snakemake.params.chunksize), snakemake.output.table)
        else:
            sys.exit(f"Request failed with code {r.status_code}: {r.text}")

    # Write query
    print("Saving query data", flush=True)
    with open(snakemake.output.query, "w") as fw:
        json.dump(data | {"query": query_txt}, fw, indent=2)

import sys
import argparse
import logging
from typing import Iterator, List
import json
import requests
import urllib.parse


SEARCH_FIELDS_URL = "https://www.ebi.ac.uk/ena/portal/api/returnFields?dataPortal=ena&result=read_run&format=json"
POST_URL = "https://www.ebi.ac.uk/ena/portal/api/search"

HEADERS = {
    "content-type" : "application/x-www-form-urlencoded",
    "connection": "keep-alive"
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


def build_query(min_date: str, max_date: str, tax_id: str, host: str, no_platform: List[str], no_lib_strat: List[str], no_lib_src: List[str], no_empty: List[str]) -> str:
    if min_date == max_date:
        date_terms = [f"collection_date={min_date}"]
    else:
        date_terms = [f"collection_date>={min_date}" if min_date else None, f"collection_date<={max_date}" if max_date else None,]
    items = [
        f"tax_tree({tax_id})" if tax_id else None,
        *date_terms,
        f'host_scientific_name="{host}"' if host else None,
        *[f'instrument_platform!="{item}"' for item in no_platform],
        *[f'library_strategy!="{item}"' for item in no_lib_strat],
        *[f'library_source!="{item}"' for item in no_lib_src],
        *[f'{field}!=""' for field in no_empty]
    ]
    return " AND ".join(item for item in items if item)


def write_iter(content: Iterator, path: str):
    with open(path, "wb") as fw:
        for chunk in content:
            fw.write(chunk)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("output")
    parser.add_argument("--start-date", type=str, required=True)
    parser.add_argument("--end-date", type=str, required=True)
    parser.add_argument("--ncbi-taxonomy", type=str, default="2697049")
    parser.add_argument("--host-scientific-name", type=str, default="Homo sapiens")
    parser.add_argument("--exclude-platform", type=str, nargs="+")
    parser.add_argument("--exclude-instrument-platform", type=str, nargs="+", default=["CAPILLARY", "DNBSEQ", "ELEMENT"])
    parser.add_argument("--exclude-library-strategy", type=str, nargs="+", default=["RNA-Seq"])
    parser.add_argument("--exclude-library-source", type=str, nargs="+", default=["TRANSCRIPTOMIC", "METAGENOMIC", "METATRANSCRIPTOMIC"])
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--chunk-size", type=int, default=1024)
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    # Init search
    logging.info("Searching ENA fields")
    fields = get_ENA_text_fields()

    # Perform search
    logging.info(f"Performing ENA search")
    data = {
        "limit": args.limit,
        "result": "read_run",
        "format": "tsv"
    }
    query_txt = build_query(
        args.start_date,
        args.end_date,
        args.ncbi_taxonomy,
        args.host_scientific_name,
        args.exclude_instrument_platform,
        args.exclude_library_strategy,
        args.exclude_library_source
    )
    data["query"] = urllib.parse.quote(query_txt, safe="")
    data["fields"] = ",".join(fields)
    logging.debug(f"Data: {data}")
    with requests.post(POST_URL, data=data, headers=HEADERS, stream=True) as r:
        if r.status_code == 200:
            logging.info("Writing ENA search results")
            write_iter(r.iter_content(chunk_size=args.chunk_size), args.output)
        else:
            sys.exit(f"Request failed with code {r.status_code}: {r.text}")

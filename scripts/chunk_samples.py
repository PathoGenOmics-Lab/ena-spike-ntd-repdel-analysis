#!/usr/bin/env python


import sys
import argparse
import csv
import itertools as it
from pathlib import Path


TEMPLATE = "{}/{}/{}/{}/{}_{}_{}"


def extract_fields(record: dict) -> list:
    return [
        record["study_accession"],
        record["sample_accession"],
        record["instrument_platform"],
        record["run_accession"],
        record["library_layout"],
        record["fastq_ftp"].count(";") + 1,
        record["library_strategy"]
    ]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--size", type=int, required=True, help="Number of items in chunk")
    parser.add_argument("output", type=Path, help="Output directory", default="chunks")
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    reader = csv.DictReader(sys.stdin, delimiter="\t")
    
    for i, records in enumerate(it.batched(reader, args.size)):
        path = args.output / f"chunk_{i}.txt"
        for record in records:
            with open(path, "w") as fw:
                fw.write(TEMPLATE.format(*extract_fields(record)) + "\n")

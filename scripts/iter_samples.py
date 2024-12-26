#!/usr/bin/env python


import sys
import argparse
import csv
import itertools as it


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
    parser.add_argument("--chunk", type=int, default=0, help="Chunk index to return, starting from 0")
    parser.add_argument("--count", action="store_true", help="Just count the number of chunks")
    args = parser.parse_args()

    reader = csv.DictReader(sys.stdin, delimiter="\t")

    if args.count:
        n = sum(1 for _ in it.batched(reader, args.size))
        sys.stdout.write(str(n))
    else:
        for i, records in enumerate(it.batched(reader, args.size)):
            if i == args.chunk:
                for record in records:
                    sys.stdout.write(TEMPLATE.format(*extract_fields(record)) + "\n")

#!/usr/bin/env python


import sys
import argparse
import csv


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
    parser.add_argument("--chunk", type=int, required=True)  # from 0
    parser.add_argument("--nrows", type=int, required=True)
    parser.add_argument("--size", type=int, required=True)
    args = parser.parse_args()

    min_row = args.size * args.chunk
    max_row = args.size * (args.chunk + 1)

    reader = csv.DictReader(sys.stdin, delimiter="\t")
    for i, record in enumerate(reader):
        if min_row <= i and i < max_row:
            sys.stdout.write(TEMPLATE.format(*extract_fields(record)) + "\n")

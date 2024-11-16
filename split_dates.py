#!/usr/bin/env python


import csv
import itertools
import argparse
import logging
from typing import List


def read_date(string: str) -> str:
    date_items = string.count("-") + 1
    if date_items == 3:
        return string
    elif date_items == 2:
        return string + "-01"
    else:
        return ""


class Group:
    def __init__(self, n: int, start: str, end: str):
        self.start = start
        self.end = end
        self.n = n
    def __add__(self, group):
        start = min(self.start, group.start)
        end = max(self.end, group.end)
        n = self.n + group.n
        return Group(n, start, end)
    def __str__(self):
        return f"Group({self.n}, {self.start}, {self.end})"
    def as_dict(self):
        return {"n": self.n, "start": self.start, "end": self.end}


def check_groups(groups: List[Group], check_size: int) -> bool:
    bound_checks = []
    for i in range(len(groups)-1):
        bound_checks.append(groups[i].end < groups[i+1].start)
    size_check = sum(group.n for group in groups) == check_size
    return all(bound_checks) and size_check


if __name__ == "__main__":

    # Parse args
    parser = argparse.ArgumentParser(description="Split ENA search results, binning approximately a certain number of records in each date bin")
    parser.add_argument("input", help="ENA search results in TSV format")
    parser.add_argument("output", help="Table containing binned dates in TSV format")
    parser.add_argument("-c", "--date-column", help="Column containing dates (assumes YYYY-MM[-DD] format)", default="collection_date")
    parser.add_argument("-n", "--n-samples", help="Approximate number of samples in each date bin", type=int, default=1000)
    parser.add_argument("--debug", help="Enable debug messages", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    logging.info(f"Reading dates from column '{args.date_column}'")
    with open(args.input) as f:
        dates = sorted(read_date(row[args.date_column]) for row in csv.DictReader(f, delimiter="\t"))
    logging.info(f"Read {len(dates)} dates")

    logging.info(f"Building groups of size around {args.n_samples}")
    group_iterator = itertools.groupby(dates)
    date, grouper = next(group_iterator)
    n = sum(1 for _ in grouper)
    groups = [Group(n, date, date)]
    for date, grouper in group_iterator:
        n = sum(1 for _ in grouper)
        group = Group(n, date, date)
        next_group = groups[-1] + group
        if next_group.n <= args.n_samples:
            logging.debug(f"Joining {groups[-1]} + {group} = {next_group}")
            groups[-1] = next_group
        else:
            logging.debug(f"Adding {group} to list")
            groups.append(group)

    logging.info(f"Writing {len(groups)} groups")
    with open(args.output, "w") as fw:
        writer = csv.DictWriter(fw, delimiter="\t", fieldnames=("start", "end", "n"))
        for group in groups:
            writer.writerow(group.as_dict())

    logging.info(f"PASS: {check_groups(groups, len(dates))}")

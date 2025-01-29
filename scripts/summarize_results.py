#!/usr/bin/env python
#SBATCH --job-name summarize
#SBATCH --mem-per-cpu 4GB
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 1
#SBATCH --qos short
#SBATCH --time 1-00:00:00
#SBATCH --output slurm-%x-%j.out


import argparse
import logging
import re
import csv
from pathlib import Path


COLUMNS = [
    "haplotype", "include_pct", "exclude_pct", "haplotype_matches", "frequency",
    "study_accession", "sample_accession", "instrument_platform", "run_accession", "library_layout", "nfastq", "library_strategy",
    "gisaid_match_column", "gisaid_pattern", "gisaid_match",
    "pango_lineage", "pango_scorpio", "pango_qc"
]
ACCNUM_P = re.compile(r"EPI_ISL_[0-9]+")
VNAME_P = re.compile(r"hCoV-19/.+/.+/[0-9]+")
RESULT_CSV_P = re.compile(r"(.+)\.inclpct_([0-9]+)\.exclpct_([0-9]+)\.csv")


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


def parse_freq(field: str) -> str:
    match field.split(","):
        case [_, _, _, _, freq, _]:
            return freq
        case _:
            return ""


def parse_results_csv(path: str) -> list:
    n = 0
    frequencies = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for record in reader:
            n += 1
            freq = parse_freq(record["S:p.Asp140_His141insValTyrTyr"]) + "|" + parse_freq(record["S:p.Val67_Ile68dup"])
            frequencies.append(freq)
    return n, ";".join(frequencies)


def parse_pangolin(path: str) -> dict:
    with open(path) as f:
        reader = csv.DictReader(f)
        # Assume 1 record
        record = next(reader)
        if record:
            logging.debug(f"Found record in '{path}'")
            return {
                "pango_lineage": record.get("lineage"),
                "pango_scorpio": record.get("scorpio_call"),
                "pango_qc": record.get("qc_status")
            }
        else:
            logging.error(f"No record in '{path}'")
            return {
                "pango_lineage": None,
                "pango_scorpio": None,
                "pango_qc": None
            }


def search_results(directory: Path) -> list:
    results = []
    for csv_path in directory.glob("*.csv"):
        result = {}
        if match := RESULT_CSV_P.match(csv_path.as_posix()):
            logging.debug(f"Matched pattern for '{csv_path}'")
            result["haplotype"], result["include_pct"], result["exclude_pct"] = match.groups()
            result["haplotype"] = Path(result["haplotype"]).name
            result["haplotype_matches"], result["frequency"] = parse_results_csv(csv_path)
            logging.debug(f"Found {result['haplotype_matches']} records in '{csv_path}'")
            results.append(result)
        else:
            logging.error(f"Could not match pattern for '{csv_path}'")
            result["haplotype"], result["include_pct"], result["exclude_pct"] = None, None, None
            result["haplotype_matches"], result["frequency"] = None, None
            results.append(result)
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("table", help="Table of read runs in TSV format")
    parser.add_argument("resdir", type=Path, help="Results directory")
    parser.add_argument("output", help="Output table in TSV format")
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    with open(args.table) as f, open(args.output, "w") as fw:
        reader = csv.DictReader(f, delimiter="\t")
        writer = csv.DictWriter(fw, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        for i, record in enumerate(reader):
            logging.debug(f"Parsing record {i}")
            fields = extract_fields(record)
            results_dir = Path(args.resdir/"repdel/filter_haplotype/{}/{}/{}/{}/{}_{}_{}".format(*fields))
            logging.debug(f"Path for record {i} is '{results_dir}'")
            # Skip if no results
            if not results_dir.is_dir():
                logging.warning(f"Directory '{results_dir}' does not exist")
                continue
            # List results on dir
            logging.debug(f"Listing results of record {i} in directory '{results_dir}'")
            results = search_results(results_dir)
            # Find GISAID identifiers
            gisaid = {"gisaid_match_column": [], "gisaid_match": [], "gisaid_pattern": []}
            for colname, value in record.items():
                for match in ACCNUM_P.findall(value):
                    gisaid["gisaid_match_column"].append(colname)
                    gisaid["gisaid_match"].append(match)
                    gisaid["gisaid_pattern"].append("accession_id")
                for match in VNAME_P.findall(value):
                    gisaid["gisaid_match_column"].append(colname)
                    gisaid["gisaid_match"].append(match)
                    gisaid["gisaid_pattern"].append("virus_name")
            logging.debug(f"Writing {len(gisaid['gisaid_pattern'])} GISAID IDs")
            gisaid = {k: ";".join(v) for k, v in gisaid.items()}
            # Find Pangolin assignment
            pangolin_path = Path(args.resdir/"pangolin/{}/{}/{}/{}/{}_{}_{}/assignment.csv".format(*fields))
            pango_reader = parse_pangolin(pangolin_path)
            # Write results
            logging.debug(f"Writing {len(results)} results of record {i} in directory '{results_dir}'")
            writer.writerows({k: v for k, v in record.items() if k in COLUMNS} | row | gisaid | pango_reader for row in results)

import csv
import argparse
import logging


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("--exclude-empty-values", type=str, nargs="+", default=["fastq_ftp"])
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    logging.info(f"Excluding empty values in {len(args.exclude_empty_values)} column(s)")
    with open(args.input) as f, open(args.output, "w") as fw:
        reader = csv.DictReader(f, delimiter="\t")
        writer = csv.DictWriter(fw, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()
        for record in reader:
            if not any(record[column] == "" for column in args.exclude_empty_values):
                writer.writerow(record)

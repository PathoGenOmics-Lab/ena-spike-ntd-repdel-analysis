import csv
import argparse
import logging


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("--subsample", type=int, default=100000)
    parser.add_argument("--omit-list", default="data/omit.csv")
    parser.add_argument("--omit-list-skip-cols", nargs="+", default="reason")
    parser.add_argument("--exclude-empty-values", nargs="+", default=["fastq_ftp"])
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    logging.info(f"Reading omitted records from {args.omit_list} excluding empty values in the {len(args.exclude_empty_values)} column(s)")
    with open(args.omit_list) as f:
        reader = csv.DictReader(f)
        omit_records = [{k: v for k, v in record.items() if k not in args.omit_list_skip_cols} for record in reader]

    logging.info("Writing subsample")
    n = 0
    with open(args.input) as f, open(args.output, "w") as fw:
        reader = csv.DictReader(f, delimiter="\t")
        writer = csv.DictWriter(fw, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()
        for record in reader:
            # Write shuffled records that pass the filter until the subsample size is reached
            if n < args.subsample and all((
                    # No empty values on any of the selected columns
                    not any(record[column] == "" for column in args.exclude_empty_values),
                    # The record is not contained in the omit list
                    not any(all(record[column] == omit[column] for column in omit.keys()) for omit in omit_records)
                )):
                writer.writerow(record)
                n += 1

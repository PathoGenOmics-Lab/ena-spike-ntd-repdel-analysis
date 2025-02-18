import sys
import csv


def count_fastq(row: dict) -> int:
    if row["fastq_ftp"]:
        return row["fastq_ftp"].count(";") + 1
    else:
        return 0


if __name__ == "__main__":
    reader = csv.DictReader(sys.stdin, delimiter="\t")
    writer = csv.DictWriter(sys.stdout, fieldnames=reader.fieldnames, delimiter="\t")
    writer.writeheader()
    paths = set()
    for row in reader:
        if count_fastq(row) in (1, 2):
            writer.writerow(row)

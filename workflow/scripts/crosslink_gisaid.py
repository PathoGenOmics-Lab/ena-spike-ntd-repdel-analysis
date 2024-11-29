import re
import csv


ACCNUM_P = re.compile(r"EPI_ISL_[0-9]+")
VNAME_P = re.compile(r"hCoV-19/.+/.+/[0-9]+")


if __name__ == "__main__":

    columns = snakemake.params.columns + ["match_column", "pattern", "match"]

    with open(snakemake.input.table) as f, open(snakemake.output.table, "w") as fw:
        reader = csv.DictReader(f, delimiter="\t")
        writer = csv.DictWriter(fw, fieldnames=columns)
        writer.writeheader()
        for record in reader:
            for colname, value in record.items():
                result = {k: v for k, v in record.items() if k in columns}
                result["match_column"] = colname
                for match in ACCNUM_P.finditer(value):
                    result["match"] = match
                    result["pattern"] = "accession_id"
                    writer.writerow(result)
                for match in VNAME_P.finditer(value):
                    result["match"] = match
                    result["pattern"] = "virus_name"
                    writer.writerow(result)

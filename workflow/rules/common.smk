SEARCH_DF_COLS = (
    "study_accession", "sample_accession", "instrument_platform",
    "run_accession", "library_layout", "fastq_ftp", "library_strategy"
)


def count_fastq(row: dict) -> int:
    return row["fastq_ftp"].count(";") + 1


def as_string(path: Path | str) -> str:
    return Path(path).as_posix()


def read_sample_paths(table: str):
    delimiter = "\t" if table.endswith(".tsv") else ","
    with open(table) as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        paths = set()
        for row in reader:
            row["nfastq"] = count_fastq(row)
            paths.add("{study_accession}/{sample_accession}/{instrument_platform}/{run_accession}/{library_layout}_{nfastq}_{library_strategy}".format(**row))
    return sorted(paths)


def read_studies(table: str):
    delimiter = "\t" if table.endswith(".tsv") else ","
    with open(table) as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        return sorted(set(row["study_accession"] for row in reader))


def read_sample_paths_from_study(table: str, study: str):
    delimiter = "\t" if table.endswith(".tsv") else ","
    with open(table) as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        paths = set()
        for row in reader:
            row["nfastq"] = count_fastq(row)
            if study == row["study_accession"]:
                paths.add("{study_accession}/{sample_accession}/{instrument_platform}/{run_accession}/{library_layout}_{nfastq}_{library_strategy}".format(**row))
    return sorted(paths)


def build_snpsift_hgvs_p_filter(wildcards):
    items = set()
    for marker_class in ("include_hgvs_p", "exclude_hgvs_p"):
        for marker in config["HAPLOTYPES"][wildcards.haplotype].get(marker_class, []):
            gene, expression = marker["gene"], marker["expression"]
            items.add(f"(ANN[*].GENE = '{gene}' & ANN[*].HGVS_P = '{expression}')")
    if len(items) > 0:
        return f" & ({' | '.join(sorted(items))})"
    else:
        return ""

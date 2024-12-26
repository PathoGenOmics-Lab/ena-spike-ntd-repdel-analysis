def count_fastq(record: dict) -> int:
    return record["fastq_ftp"].count(";") + 1


def iter_records(cursor):
    colnames = [field[0] for field in cursor.description]
    for row_values in cursor:
        yield {col: value for col, value in zip(colnames, row_values)}


def read_sample_paths(database: str):
    with sqlite3.connect(database, timeout=30) as conn:
        cursor = conn.execute("SELECT study_accession, sample_accession, instrument_platform, run_accession, library_layout, library_strategy, fastq_ftp FROM ENARecords")
        if cursor is None:
            return []
        paths = set()
        for record in iter_records(cursor):
            print(record)
            record["nfastq"] = count_fastq(record)
            paths.add("{study_accession}/{sample_accession}/{instrument_platform}/{run_accession}/{library_layout}_{nfastq}_{library_strategy}".format(**record))
    return sorted(paths)


def read_sample_paths_from_study(database: str, study: str):
    with sqlite3.connect(database, timeout=30) as conn:
        cursor = conn.execute("SELECT study_accession, sample_accession, instrument_platform, run_accession, library_layout, library_strategy, fastq_ftp FROM ENARecords")
        if cursor is None:
            return []
        paths = set()
        for record in iter_records(cursor):
            record["nfastq"] = count_fastq(record)
            if study == record["study_accession"]:
                paths.add("{sample_accession}/{instrument_platform}/{run_accession}/{library_layout}_{nfastq}_{library_strategy}".format(**record))
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


rule cat_csv:
    threads: 1
    input: ["placeholder"]
    output: "placeholder"
    resources:
        runtime = "15m",
        mem_mb = 500
    run:
        import csv
        import logging
        logging.basicConfig(
            level=logging.INFO,
            format=config["PY_LOG_FMT"],
            filename=log[0]
        )
        # Read header from first input
        with open(input[0]) as f:
            header = csv.DictReader(f).fieldnames
        # Write records
        n = 0
        with open(output[0], "w") as fw:
            writer = csv.DictWriter(fw, fieldnames=header)
            writer.writeheader()
            for table_path in input:
                with open(table_path) as f:
                    reader = csv.DictReader(f)
                    for record in reader:
                        writer.writerow(record)
                        n += 1
        logging.info(f"Wrote {n} records")

rule consensus_merge:
    group: "pangolin_{study}"
    input: expand(OUTPUT/"variants/consensus/{path}/sample.fasta", path=SAMPLE_PATHS)
    output: temp(OUTPUT/"pangolin/sequences.fasta")
    resources:
        runtime = "15m",
        mem_gb = 2
    shell: "cat {input:q} > {output:q}"


rule pangolin_assignment:
    threads: 32
    group: "pangolin_{study}"
    conda: "../envs/lineages.yaml"
    shadow: "minimal"
    input:
        fasta = OUTPUT/"pangolin/sequences.fasta"
    output:
        table = OUTPUT/"pangolin/pangolin.csv"
    resources:
        mem_gb = 8,
        max_cpu_per_node = lambda wc, threads: threads
    log: OUTPUT/"logs/pangolin/pangolin_assignment.txt"
    shell: "pangolin {input.fasta:q} --outfile {output.table:q} --threads {threads} >{log:q} 2>&1"


rule filter_pangolin:
    input:
        table = OUTPUT/"pangolin/pangolin.csv"
    params:
        qc_status = ["pass"],
        scorpio_call = ["Omicron (BA.1-like)"],
        lineage_base = ["BA.1"]
    output:
        table = temp(OUTPUT/"pangolin/pangolin.filtered.csv")
    log: OUTPUT/"logs/pangolin/filter_pangolin.txt"
    resources:
        runtime = "20m",
        mem_gb = 2
    run:
        import logging
        import csv
        logging.basicConfig(
            level=logging.INFO,
            format=config["PY_LOG_FMT"],
            filename=log[0]
        )
        n = 0
        ntotal = 0
        with open(input.table) as f, open(output.table, "w") as fw:
            reader = csv.DictReader(f)
            writer = csv.DictWriter(fw, fieldnames=reader.fieldnames)
            writer.writeheader()
            for row in reader:
                if row["qc_status"] in params.qc_status and any((
                    row["scorpio_call"] in params.scorpio_call,
                    any(row["lineage"].startswith(lineage) for lineage in params.lineage_base)
                )):
                    writer.writerow(row)
                    n += 1
                ntotal += 1
        logging.info(f"Wrote {n} of {ntotal} records with qc_status in {params.qc_status}, and scorpio_call in {params.scorpio_call} or lineage base in {params.lineage_base}")


checkpoint select_samples_after_processing:
    input:
        search_table = config["SEARCH_TABLE"],
        coverage_table = OUTPUT/"variants/coverage.filtered.csv",
        pangolin_table = OUTPUT/"pangolin/pangolin.filtered.csv"
    params:
        consensus_pattern = r"Consensus_([A-Z0-9]+)__([A-Z0-9]+)__([A-Z_]+)__([A-Z0-9]+)__([A-Z]+)__([A-Z]+)_threshold_[0-9\.]+_quality_[0-9\.]+",
        coverage_sample_pattern = r"([A-Z0-9]+)__([A-Z0-9]+)__([A-Z_]+)__([A-Z0-9]+)__([A-Z]+)__([A-Z]+)"
    output:
        search_table = OUTPUT/"ena/search.filtered.afterproc.tsv"
    resources:
        mem_mb = 16000,
        runtime = "30m"
    log: OUTPUT/"logs/pangolin/select_samples_after_processing.txt"
    script: "../scripts/select_samples_after_processing.py"

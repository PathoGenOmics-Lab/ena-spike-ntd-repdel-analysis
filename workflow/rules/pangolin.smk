rule pangolin_assignment:
    threads: 1
    conda: "../envs/lineages.yaml"
    shadow: "minimal"
    input:
        fasta = OUTPUT/"variants/consensus/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.fasta"
    output:
        table = OUTPUT/"pangolin/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/assignment.csv"
    resources:
        mem_mb = 4000
    log: OUTPUT/"logs/pangolin/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/assignment.txt"
    shell: "pangolin {input.fasta:q} --outfile {output.table:q} --threads {threads} >{log:q} 2>&1"


rule filter_pangolin:
    input:
        table = OUTPUT/"pangolin/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/assignment.csv"
    params:
        qc_status = ["pass"],
        scorpio_call = ["Omicron (BA.1-like)"],
        lineage_base = ["BA.1"]
    output:
        table = temp(OUTPUT/"pangolin/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/assignment.filtered.csv")
    log: OUTPUT/"logs/pangolin/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/filter.txt"
    resources:
        runtime = "10m",
        mem_mb = 250
    run:
        import logging
        import csv
        logging.basicConfig(
            level=logging.INFO,
            format=config["PY_LOG_FMT"],
            filename=log[0]
        )
        passed = False
        with open(input.table) as f, open(output.table, "w") as fw:
            reader = csv.DictReader(f)
            writer = csv.DictWriter(fw, fieldnames=reader.fieldnames)
            writer.writeheader()
            row = next(reader)
            if row["qc_status"] in params.qc_status and any((
                row["scorpio_call"] in params.scorpio_call,
                any(row["lineage"].startswith(lineage) for lineage in params.lineage_base)
            )):
                writer.writerow(row)
                passed = True
        logging.info(f"Record pass={passed}")


rule pangolin_merge:
    shadow: "shallow"
    input: expand(OUTPUT/"pangolin/{path}/assignment.filtered.csv", path=SAMPLE_PATHS)
    output: OUTPUT/"variants/pangolin.filtered.csv"
    resources:
        runtime = "10m",
        mem_mb = 500
    log: OUTPUT/"logs/pangolin/pangolin_merge.txt"
    shell: "head -n 1 {input[0]:q} >{output:q} 2>{log:q} && tail -n +2 -q {input:q} >>{output:q} 2>>{log:q}"



rule select_samples_after_processing:
    input:
        search_table = config["SEARCH_TABLE"],
        coverage_table = OUTPUT/"variants/coverage.filtered.csv",
        pangolin_table = OUTPUT/"pangolin/pangolin.filtered.csv"
    params:
        consensus_pattern = r"Consensus_([A-Z0-9]+)__([A-Z0-9]+)__([A-Z_]+)__([A-Z0-9]+)__([A-Z]+)__([A-Z]+)_threshold_[0-9\.]+_quality_[0-9\.]+",
        coverage_sample_pattern = r"([A-Z0-9]+)__([A-Z0-9]+)__([A-Z_]+)__([A-Z0-9]+)__([A-Z]+)__([A-Z]+)"
    output:
        search_table = config["FILTERED_TABLE"]
    resources:
        mem_mb = 16000,
        runtime = "30m"
    log: OUTPUT/"logs/pangolin/select_samples_after_processing.txt"
    script: "../scripts/select_samples_after_processing.py"

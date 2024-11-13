rule consensus_merge:
    group: "pangolin_{study}"
    input: lambda w: build_search_targets_filtering(w, f"output/variants/consensus/{w.study}/{{}}/{{}}/{{}}/{{}}_{{}}_{{}}/sample.fasta", ("sample_accession", "instrument_platform", "run_accession", "library_layout", "fastq_ftp", "library_strategy"), study_accession=w.study)
    output: temp("output/pangolin/consensus_merge/{study}/sequences.fasta")
    resources:
        runtime = "15m",
        mem_gb = 2
    shell: "cat {input:q} > {output:q}"


rule pangolin_assignment:
    threads: 8
    group: "pangolin_{study}"
    conda: "../envs/lineages.yaml"
    shadow: "minimal"
    input:
        fasta = "output/pangolin/consensus_merge/{study}/sequences.fasta"
    output:
        table = temp("output/pangolin/pangolin_assignment/{study}/pangolin.csv")
    resources:
        mem_gb = 8
    log: "output/logs/pangolin/pangolin_assignment/{study}.txt"
    shell: "pangolin {input.fasta:q} --outfile {output.table:q} --threads {threads} >{log:q} 2>&1"


rule pangolin_assignment_merge:
    input: lambda w: build_search_targets(w, "output/pangolin/pangolin_assignment/{}/pangolin.csv", ("study_accession",))
    output: "output/pangolin/pangolin.csv"
    resources:
        runtime = "15m",
        mem_gb = 2
    shell: "head -n 1 {input[0]:q} >{output:q} && tail -n +2 -q {input:q} >>{output:q}"


rule filter_pangolin:
    input:
        table = "output/pangolin/pangolin.csv"
    params:
        qc_status = ["pass"],
        scorpio_call = ["Omicron (BA.1-like)"],
        lineage_base = ["BA.1"]
    output:
        table = "output/pangolin/pangolin.filtered.csv"
    log: "output/logs/pangolin/filter_pangolin.txt"
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


checkpoint filter_search_ena_with_pangolin:
    input:
        search_table = "output/ena/search.filtered.tsv",
        pangolin_table = "output/pangolin/pangolin.filtered.csv"
    params:
        # Consensus_{study}__{sample}__{platform}__{run}__{layout}__{strategy}_threshold_{threshold}_quality_{quality}
        pattern = r"Consensus_([A-Z0-9]+)__([A-Z0-9]+)__([A-Z_]+)__([A-Z0-9]+)__([A-Z]+)__([A-Z]+)_threshold_[0-9\.]+_quality_[0-9\.]+"
    output:
        search_table = "output/ena/search.filtered.pangolin.tsv"
    resources:
        mem_mb = 16000,
        runtime = "30m"
    log: "output/logs/pangolin/filter_search_ena_with_pangolin.txt"
    run:
        import logging
        import re
        import pandas as pd
        logging.basicConfig(
            level=logging.INFO,
            format=config["PY_LOG_FMT"],
            filename=log[0]
        )
        columns = [
            "study_accession", "sample_accession", "instrument_platform",
            "run_accession", "library_layout", "library_strategy"
        ]
        id_pattern = re.compile(params.pattern)
        # Read pangolin assignment
        logging.info("Reading pangolin assignment")
        pangolin = pd.read_csv(input.pangolin_table)
        logging.info(f"Read {len(pangolin)} records")
        # Extract sample metadata from record ID
        logging.info("Extracting metadata from taxons")
        if len(pangolin) != 0:
            pangolin[columns] = pangolin["taxon"].apply(lambda record_id: id_pattern.match(record_id).groups()).tolist()
            # Keep only sample metadata
            pangolin = pangolin.drop([col for col in columns if col not in columns], axis="columns")
        # Read filtered search results
        logging.info("Reading filtered search results")
        search = pd.read_csv(input.search_table, sep="\t")
        # Keep only search results in the pangolin table and write
        logging.info("Filtering and writing")
        search.merge(
            pangolin,
            how="right",
            on=columns
        ).to_csv(output.search_table, sep="\t", index=False)

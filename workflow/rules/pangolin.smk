rule consensus_merge:
    input: lambda w: build_search_targets_filtering(w, f"output/variants/consensus/{w.study}/{{}}/{{}}/{{}}/{{}}_{{}}_{{}}/sample.fasta", ("sample_accession", "instrument_platform", "run_accession", "library_layout", "fastq_ftp", "library_strategy"), study_accession=w.study)
    output: temp("output/pangolin/consensus_merge/{study}/sequences.fasta")
    resources:
        runtime = "15m",
        mem_gb = 2
    shell: "cat {input:q} > {output:q}"


rule pangolin_assignment:
    threads: 8
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


checkpoint filter_pangolin:
    input:
        table = "output/pangolin/pangolin.csv"
    params:
        qc_status = ["pass"],
        scorpio_call = ["Omicron (BA.1-like)"]
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
                if row["qc_status"] in params.qc_status and row["scorpio_call"] in params.scorpio_call:
                    writer.writerow(row)
                    n += 1
                ntotal += 1
        logging.info(f"Wrote {n} of {ntotal} records with qc_status in {params.qc_status} and scorpio_call in {params.scorpio_call}")

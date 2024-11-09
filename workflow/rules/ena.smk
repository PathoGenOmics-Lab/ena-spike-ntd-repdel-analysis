rule search_ena:
    params:
        start_date = "2021-11-01",       # outbreak.info approx BA.1 start date
        end_date = "2022-08-01",         # outbreak.info approx BA.1 end date
        limit = config["SEARCH_LIMIT"],  # 0 means no record limit
        taxonomy = "2697049",            # SARS-CoV-2
        host_scientific_name = "Homo sapiens",
        chunksize = 1024
    output:
        table = "output/ena/search.tsv",
        query = "output/ena/search.json"
    log: "output/logs/ena/search_ena.txt"
    script: "../scripts/search_ena.py"


checkpoint filter_search_ena:
    input:
        table = "output/ena/search.tsv"
    params:
        omit_platform = ["CAPILLARY", "DNBSEQ", "ELEMENT"],
        omit_library_strategy = ["RNA-Seq"]
    output:
        table = "output/ena/search.filtered.tsv",
    resources:
        mem_mb = 12000,
        runtime = "15m"
    run:
        import pandas as pd
        df = pd.read_csv(input.table, sep="\t")
        df[
            ~df["instrument_platform"].isin(params.omit_platform) & \
            ~df["library_strategy"].isin(params.omit_library_strategy) & \
            df["fastq_ftp"].str.count(";").isin([0, 1])
        ].to_csv(output.table, sep="\t", index=False)


rule summarize_ena_search:
    conda: "../envs/rdata.yaml"
    input:
        table = "output/ena/search.filtered.tsv"
    params:
        count_bins = 9,
        plot_width_in = 25
    resources:
        mem_mb = 12000,
        runtime = "15m"
    output:
        plot_pdf = "output/ena/summarize_ena_search/summary.pdf",
        plot_png = "output/ena/summarize_ena_search/summary.png",
        country_timeline_table = "output/ena/summarize_ena_search/summary_country_timeline.csv",
        tech_timeline_table = "output/ena/summarize_ena_search/summary_tech_timeline.csv",
        seqres_timeline_table = "output/ena/summarize_ena_search/summary_seqres_timeline.csv",
        studies_table = "output/ena/summarize_ena_search/studies.csv"
    log: "output/logs/ena/summarize_ena_search.txt"
    script: "../scripts/summarize_ena_search.R"


use rule summarize_ena_search as summarize_ena_search_pangolin with:
    input:
        table = "output/ena/search.filtered.pangolin.tsv"
    output:
        plot_pdf = "output/ena/summarize_ena_search_pangolin/summary.pdf",
        plot_png = "output/ena/summarize_ena_search_pangolin/summary.png",
        country_timeline_table = "output/ena/summarize_ena_search_pangolin/summary_country_timeline.csv",
        tech_timeline_table = "output/ena/summarize_ena_search_pangolin/summary_tech_timeline.csv",
        seqres_timeline_table = "output/ena/summarize_ena_search_pangolin/summary_seqres_timeline.csv",
        studies_table = "output/ena/summarize_ena_search_pangolin/studies.csv"
    log: "output/logs/ena/summarize_ena_search_pangolin.txt"


rule split_ena_search_results:
    group: "group_ena_{study}"
    input:
        table = "output/ena/search.filtered.tsv"
    output:
        table = "output/ena/search/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/runs.csv"
    resources:
        runtime = "10m",
        mem_gb = 2
    log: "output/logs/ena/split_ena_search_results/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}.txt"
    run:
        import logging
        import csv
        logging.basicConfig(
            level=logging.INFO,
            format=config["PY_LOG_FMT"],
            filename=log[0]
        )
        n = 0
        with open(input.table) as f, open(output.table, "w") as fw:
            reader = csv.DictReader(f, delimiter="\t")
            writer = csv.DictWriter(fw, fieldnames=reader.fieldnames)
            writer.writeheader()
            for row in reader:
                if all((
                    row["run_accession"] == wildcards.run,
                    row["sample_accession"] == wildcards.sample,
                    row["study_accession"] == wildcards.study,
                    row["instrument_platform"] == wildcards.platform,
                    row["library_layout"] == wildcards.layout,
                    row["library_strategy"] == wildcards.strategy
                )):
                    writer.writerow(row)
                    n += 1
        logging.info(f"Wrote {n} records with study={wildcards.study}, sample={wildcards.sample}, platform={wildcards.platform}, run={wildcards.run}, layout={wildcards.layout} and strategy={wildcards.strategy}")


rule download_ena_one_fastq:
    group: "group_ena_{study}"
    input:
        table = "output/ena/search/{study}/{sample}/{platform}/{run}/{layout}_1_{strategy}/runs.csv"
    params:
        retries = 5,
        backoff_factor = 1,
        backoff_jitter = 1
    output:
        fastq = temp("output/ena/downloads/fastq/{study}/{sample}/{platform}/{run}/{layout}_1_{strategy}/sample.fastq.gz")
    log: "output/logs/ena/download_ena/{study}/{sample}/{platform}/{run}/{layout}_1_{strategy}.txt"
    resources:
        ena_api_calls_per_second = 1,
        runtime = "30m"
    retries: 2
    script: "../scripts/download_ena_one_fastq.py"


rule download_ena_two_fastq:
    group: "group_ena_{study}"
    input:
        table = "output/ena/search/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/runs.csv"
    params:
        retries = 5,
        backoff_factor = 1,
        backoff_jitter = 1
    output:
        fastq_1 = temp("output/ena/downloads/fastq/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/sample.R1.fastq.gz"),
        fastq_2 = temp("output/ena/downloads/fastq/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/sample.R2.fastq.gz")
    resources:
        ena_api_calls_per_second = 1,
        runtime = "30m"
    retries: 2
    log: "output/logs/ena/download_ena/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}.txt"
    script: "../scripts/download_ena_two_fastq.py"

rule summarize_ena_search:
    conda: "../envs/rdata.yaml"
    input:
        table = config["SEARCH_TABLE"]
    params:
        count_bins = 9,
        plot_width_in = 25
    resources:
        mem_mb = 12000,
        runtime = "15m"
    output:
        plot_pdf = OUTPUT/"ena/summarize_ena_search/summary.pdf",
        plot_png = OUTPUT/"ena/summarize_ena_search/summary.png",
        country_timeline_table = OUTPUT/"ena/summarize_ena_search/summary_country_timeline.csv",
        tech_timeline_table = OUTPUT/"ena/summarize_ena_search/summary_tech_timeline.csv",
        seqres_timeline_table = OUTPUT/"ena/summarize_ena_search/summary_seqres_timeline.csv",
        studies_table = OUTPUT/"ena/summarize_ena_search/studies.csv",
        tech_table = OUTPUT/"ena/summarize_ena_search/technologies.csv"
    log: OUTPUT/"logs/ena/summarize_ena_search.txt"
    script: "../scripts/summarize_ena_search.R"


use rule summarize_ena_search as summarize_ena_search_after_processing with:
    input:
        table = config["FILTERED_TABLE"]
    output:
        plot_pdf = OUTPUT/"ena/summarize_ena_search_after_processing/summary.pdf",
        plot_png = OUTPUT/"ena/summarize_ena_search_after_processing/summary.png",
        country_timeline_table = OUTPUT/"ena/summarize_ena_search_after_processing/summary_country_timeline.csv",
        tech_timeline_table = OUTPUT/"ena/summarize_ena_search_after_processing/summary_tech_timeline.csv",
        seqres_timeline_table = OUTPUT/"ena/summarize_ena_search_after_processing/summary_seqres_timeline.csv",
        studies_table = OUTPUT/"ena/summarize_ena_search_after_processing/studies.csv",
        tech_table = OUTPUT/"ena/summarize_ena_search_after_processing/technologies.csv"
    log: OUTPUT/"logs/ena/summarize_ena_search_after_processing.txt"


rule split_ena_search_results:
    group: "download"
    input:
        table = config["SEARCH_TABLE"]
    output:
        table = temp(OUTPUT/"ena/search/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/runs.csv")
    resources:
        runtime = lambda wc, attempt: 15 * attempt,
        mem_mb = lambda wc, attempt: 4000 * attempt
    retries: 2
    log: OUTPUT/"logs/ena/split_ena_search_results/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}.txt"
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
    group: "download"
    input:
        table = OUTPUT/"ena/search/{study}/{sample}/{platform}/{run}/{layout}_1_{strategy}/runs.csv"
    params:
        retries = 5,
        backoff_factor = 1,
        backoff_jitter = 1,
        chunk_mb = 10
    output:
        fastq = temp(OUTPUT/"ena/downloads/fastq/{study}/{sample}/{platform}/{run}/{layout}_1_{strategy}/sample.fastq.gz")
    log: OUTPUT/"logs/ena/download_ena/{study}/{sample}/{platform}/{run}/{layout}_1_{strategy}.txt"
    resources:
        ena_api_calls_per_second = 1,
        runtime = lambda wc, attempt: 30 * attempt,
        mem_mb = lambda wc, attempt: 2000 * attempt
    retries: 2
    script: "../scripts/download_ena_one_fastq.py"


rule download_ena_two_fastq:
    group: "download"
    input:
        table = OUTPUT/"ena/search/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/runs.csv"
    params:
        retries = 5,
        backoff_factor = 1,
        backoff_jitter = 1,
        chunk_mb = 10
    output:
        fastq_1 = temp(OUTPUT/"ena/downloads/fastq/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/sample.R1.fastq.gz"),
        fastq_2 = temp(OUTPUT/"ena/downloads/fastq/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/sample.R2.fastq.gz")
    resources:
        ena_api_calls_per_second = 1,
        runtime = lambda wc, attempt: 30 * attempt,
        mem_mb = lambda wc, attempt: 2000 * attempt
    retries: 2
    log: OUTPUT/"logs/ena/download_ena/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}.txt"
    script: "../scripts/download_ena_two_fastq.py"


rule crosslink_gisaid:
    threads: 1
    input:
        table = config["SEARCH_TABLE"]
    params:
        columns = ["study_accession", "sample_accession", "instrument_platform", "run_accession", "library_layout", "library_strategy"]
    output:
        table = OUTPUT/"ena/crosslink_gisaid.csv"
    resources:
        runtime = "60m",
        mem_mb = 4000
    log: OUTPUT/"logs/ena/crosslink_gisaid.txt"
    script: "../scripts/crosslink_gisaid.py"

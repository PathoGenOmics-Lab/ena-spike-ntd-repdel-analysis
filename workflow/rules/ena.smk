checkpoint search_ena:
    conda: "../envs/pydata.yaml"
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


rule filter_search_ena:
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
            ~df["library_strategy"].isin(params.omit_library_strategy)
        ].to_csv(output.table, sep="\t")


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
        plot_pdf = "output/ena/report/search/summary.pdf",
        plot_png = "output/ena/report/search/summary.png",
        country_timeline_table = "output/ena/report/search/summary_country_timeline.csv",
        tech_timeline_table = "output/ena/report/search/summary_tech_timeline.csv",
        seqres_timeline_table = "output/ena/report/search/summary_seqres_timeline.csv"
    log: "output/logs/ena/summarize_ena_search.txt"
    script: "../scripts/summarize_ena_search.R"


rule split_ena_search_results:
    group: "download"
    input:
        table = "output/ena/search.filtered.tsv"
    output:
        table = "output/ena/search/{study}/{sample}/{platform}/{run}/{layout}_{strategy}/runs.csv"
    resources:
        runtime = "10m",
        mem_mb = 12000
    run:
        import pandas as pd
        df = pd.read_csv(input.table, sep="\t")
        df[
            (df["instrument_platform"] == wildcards.platform) & \
            (df["run_accession"] == wildcards.run) & \
            (df["sample_accession"] == wildcards.sample) & \
            (df["study_accession"] == wildcards.study) & \
            (df["library_layout"] == wildcards.layout)
        ].to_csv(output.table)


rule download_ena:
    group: "download"
    conda: "../envs/pydata.yaml"
    input:
        table = "output/ena/search/{study}/{sample}/{platform}/{run}/{layout}_{strategy}/runs.csv"
    params:
        sleep = 1
    output:
        folder = directory("output/ena/downloads/fastq/{study}/{sample}/{platform}/{run}/{layout}_{strategy}")
    log: "output/logs/ena/download_ena/{study}/{sample}/{platform}/{run}/{layout}_{strategy}.txt"
    script: "../scripts/download_ena.py"

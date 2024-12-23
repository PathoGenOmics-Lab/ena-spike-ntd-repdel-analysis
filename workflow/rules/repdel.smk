rule filter_haplotype:
    conda: "../envs/pydata.yaml"
    input:
        variants = OUTPUT/"variants/snpsift_extract_variants/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/{haplotype}.tsv",
        pangolin = OUTPUT/"pangolin/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/assignment.filtered.csv"
    params:
        columns = ["CHROM", "REF", "POS", "ALT", "DP", "ALT_DP", "ALT_RV", "ALT_FREQ", "ALT_QUAL", "GENE", "HGVS_P"],
        markers = lambda w: config["HAPLOTYPES"][w.haplotype]
    output:
        # inclpct: int - minimum frequency threshold (%) to pass an "included" marker
        # exclpct: int - maximum frequency threshold (%) to pass an "excluded" marker
        table = OUTPUT/"repdel/filter_haplotype/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/{haplotype}.inclpct_{inclpct}.exclpct_{exclpct}.csv"
    log: OUTPUT/"logs/repdel/filter_haplotype/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/{haplotype}.inclpct_{inclpct}.exclpct_{exclpct}.txt"
    script: "../scripts/filter_haplotype.py"


use rule cat_csv as merge_haplotypes with:
    input: expand(OUTPUT/"repdel/filter_haplotype/{path}/{{haplotype}}.inclpct_{{inclpct}}.exclpct_{{exclpct}}.csv", path=SAMPLE_PATHS)
    output: OUTPUT/"repdel/merge_haplotypes/{haplotype}.inclpct_{inclpct}.exclpct_{exclpct}.csv"
    log: OUTPUT/"logs/repdel/merge_haplotypes/{haplotype}.inclpct_{inclpct}.exclpct_{exclpct}.txt"


rule report_region:
    conda: "../envs/rdata.yaml"
    input:
        bam = OUTPUT/"mapping/sorted_bam/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.sorted.bam",
        reference = "data/snpEff/data/{}/sequences.fa".format(config["REFERENCE"])
    params:
        plot_width_per_position_in = 0.1,
        plot_height_in = 10,
        positions = list(range(config["REPORT"]["START"], config["REPORT"]["END"]+1)),
        max_depth = 2e9,
        min_base_quality = 0,
        min_mapq = 0,
        min_nucleotide_depth = 0,
        min_minor_allele_depth = 0,
        distinguish_strands = False,
        distinguish_nucleotides = True,
        ignore_query_Ns = True,
        include_deletions = True,
        include_insertions = False
    output:
        plot = OUTPUT/"repdel/report_region/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.png",
        plot_data = OUTPUT/"repdel/report_region/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.csv"
    log: OUTPUT/"logs/repdel/report_region/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}.txt"
    script: "../scripts/report_region.R"

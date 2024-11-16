checkpoint filter_haplotype:
    conda: "../envs/pydata.yaml"
    input:
        lambda w: build_afterproc_targets(w, fOUTPUT/"variants/snpsift_extract_variants/{{}}/{{}}/{{}}/{{}}/{{}}_{{}}_{{}}/{w.haplotype}.tsv")
    params:
        columns = ["CHROM", "REF", "POS", "ALT", "DP", "ALT_DP", "ALT_RV", "ALT_FREQ", "ALT_QUAL", "GENE", "HGVS_P"],
        markers = lambda w: config["HAPLOTYPES"][w.haplotype]
    output:
        # inclpct: int - minimum frequency threshold (%) to pass an "included" marker
        # exclpct: int - maximum frequency threshold (%) to pass an "excluded" marker
        table = OUTPUT/"repdel/filter_haplotype/{haplotype}.inclpct_{inclpct}.exclpct_{exclpct}.csv"
    log: OUTPUT/"logs/repdel/filter_haplotype/{haplotype}.inclpct_{inclpct}.exclpct_{exclpct}.txt"
    script: "../scripts/filter_haplotype.py"


rule report_region:
    conda: "../envs/rdata.yaml"
    input:
        bam = OUTPUT/"mapping/sorted_bam/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.sorted.bam"
    params:
        plot_width_per_position_in = 0.5,
        plot_height_in = 10,
        positions = list(range(config["COVERAGE_FILTER"]["START"], config["COVERAGE_FILTER"]["END"]+1)),
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

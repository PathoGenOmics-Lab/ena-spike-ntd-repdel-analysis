rule filter_haplotype:
    conda: "../envs/pydata.yaml"
    input:
        lambda w: build_pangolin_targets(w, f"output/variants/snpsift_extract_variants/{{study}}/{{sample}}/{{platform}}/{{run}}/{{layout}}_{{nfastq}}_{{strategy}}/{w.haplotype}.tsv")
    params:
        columns = ["CHROM", "REF", "POS", "ALT", "DP", "ALT_DP", "ALT_RV", "ALT_FREQ", "ALT_QUAL", "GENE", "HGVS_P"],
        markers = lambda w: config["HAPLOTYPES"][w.haplotype]
    output:
        # inclpct: int - minimum frequency threshold (%) to pass an "included" marker
        # exclpct: int - maximum frequency threshold (%) to pass an "excluded" marker
        table = "output/repdel/filter_haplotype/{haplotype}.inclpct_{inclpct}.exclpct_{exclpct}.csv"
    log: "output/logs/repdel/filter_haplotype/{haplotype}.inclpct_{inclpct}.exclpct_{exclpct}.txt"
    script: "../scripts/filter_haplotype.py"


rule report_region:
    conda: "../envs/rdata.yaml"
    input:
        tables = lambda w: build_pangolin_targets(w, "output/variants/variant_calling/{}/{}/{}/{}/{}_{}_{}/sample.tsv")
    params:
        plot_width_per_position_in = 0.5,
        plot_height_in = 10,
        positions = [list(range(21765, 21771)) + list(range(21991, 21994)) + list(range(21987, 21996)), list(range(21990, 21999))],
        filter_pass = True
    output:
        plot = "output/repdel/report_region/region.png",
        plot_data = "output/repdel/report_region/region.csv"
    log: "output/logs/repdel/report_region.txt"
    script: "../scripts/report_region.R"

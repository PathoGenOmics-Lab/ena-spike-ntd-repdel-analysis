rule report_region:
    conda: "../envs/rdata.yaml"
    input:
        tables = lambda w: build_pangolin_targets(w, "output/variants/variant_calling/{}/{}/{}/{}/{}_{}_{}/sample.tsv")
    params:
        plot_width_per_position_in = 0.5,
        plot_height_in = 15,
        positions = [list(range(21765, 21771)) + list(range(21991, 21994)) + list(range(21987, 21996)), list(range(21990, 21999))],
        filter_pass = True
    output:
        plot = "output/repdel/report_region/region.png",
        plot_data = "output/repdel/report_region/region.csv"
    log: "output/logs/repdel/report_region.txt"
    script: "../scripts/report_region.R"

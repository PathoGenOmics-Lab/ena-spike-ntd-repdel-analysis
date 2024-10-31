rule variant_calling:
    conda: "../envs/reads.yaml"
    shadow: "minimal"
    params:
        max_depth = 0,
        min_quality = 0,
        ivar_quality = 20,
        ivar_freq = 0.05,
        ivar_depth = 30
    input:
        reference = "output/reference/sequence.fasta",
        features = "output/reference/features.filtered.gff"
    shell:
        "samtools mpileup -aa -x -A -d {params.max_depth} -B -Q {params.min_quality} -f {input.reference:q} {input.bam:q} | "
        "ivar variants -p result -q {params.ivar_quality} -t {params.ivar_freq} -m {params.ivar_depth} -g {input.gff} -r {input.reference:q} -g {input.features:q} &&"
        "mv result.tsv {output.tsv:q}"

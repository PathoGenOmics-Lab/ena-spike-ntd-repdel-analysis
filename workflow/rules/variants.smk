rule pileup:
    conda: "../envs/reads.yaml"
    input:
        reference = "output/reference/sequence.fasta",
        bam = "output/mapping/sorted_bam/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.sorted.bam"
    params:
        max_depth = 0,   # 0 means unrestricted
        min_quality = 0  # filtered later with iVar
    output:
        pileup = "output/variants/pileup/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.pileup"
    log: "output/logs/variants/pileup/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}.txt"
    shell: "samtools mpileup -aa -x -A -d {params.max_depth} -B -Q {params.min_quality} -f {input.reference:q} {input.bam:q} >{output.pileup:q} 2>{log}"


rule consensus:
    conda: "../envs/reads.yaml"
    shadow: "minimal"
    input:
        pileup = "output/variants/pileup/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.pileup"
    params:
        min_quality = 20,
        min_frequency = 0,
        min_ins_frequency = 0.9,
        min_depth = 30,
        char_under_min_depth = "N"
    output:
        fasta = "output/variants/consensus/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.fasta",
        quality = "output/variants/consensus/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.qual.txt"
    log: "output/logs/variants/consensus/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}.txt"
    shell: "ivar consensus -p result -q {params.min_quality} -t {params.min_frequency} -m {params.min_depth} -c {params.min_ins_frequency} -n {params.char_under_min_depth} <{input.pileup:q} 2>{log} && mv result.fa {output.fasta:q} && mv result.qual.txt {output.quality:q}"


rule variant_calling:
    conda: "../envs/reads.yaml"
    shadow: "minimal"
    params:
        min_quality = 20,
        min_frequency = 0.05,
        min_depth = 30
    input:
        reference = "output/reference/sequence.fasta",
        gff = "output/reference/features.filtered.gff3",
        pileup = "output/variants/pileup/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.pileup"
    output:
        tsv = "output/variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.tsv"
    log: "output/logs/variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}.txt"
    shell: "ivar variants -p result -q {params.min_quality} -t {params.min_frequency} -m {params.min_depth} -g {input.gff} -r {input.reference:q} <{input.pileup:q} 2>{log} && mv result.tsv {output.tsv:q}"

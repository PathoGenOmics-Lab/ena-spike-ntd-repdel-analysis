rule pileup:
    group: "sample"
    conda: "../envs/reads.yaml"
    input:
        reference = "data/snpEff/data/{}/sequences.fa".format(config["REFERENCE"]),
        bam = OUTPUT/"mapping/sorted_bam/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.sorted.bam"
    params:
        max_depth = 0,   # 0 means unrestricted
        min_quality = 0  # filtered later with iVar
    output:
        pileup = temp(OUTPUT/"variants/pileup/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.pileup")
    resources:
        runtime = lambda wc, attempt: 20 * attempt,
        mem_mb = lambda wc, attempt: 2000 * attempt
    retries: 2
    log: OUTPUT/"logs/variants/pileup/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}.txt"
    shell: "samtools mpileup -aa -x -A -d {params.max_depth} -B -Q {params.min_quality} -f {input.reference:q} {input.bam:q} >{output.pileup:q} 2>{log:q}"


rule consensus:
    group: "sample"
    conda: "../envs/reads.yaml"
    shadow: "minimal"
    input:
        pileup = OUTPUT/"variants/pileup/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.pileup"
    params:
        min_quality = 20,
        min_frequency = 0,
        min_ins_frequency = 0.9,
        min_depth = 30,
        char_under_min_depth = "N",
        prefix = "{study}__{sample}__{platform}__{run}__{layout}__{strategy}"
    output:
        fasta = OUTPUT/"variants/consensus/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.fasta",
        quality = OUTPUT/"variants/consensus/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.qual.txt"
    log: OUTPUT/"logs/variants/consensus/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}.txt"
    resources:
        runtime = lambda wc, attempt: 30 * attempt,
        mem_mb = lambda wc, attempt: 4000 * attempt
    retries: 2
    shell:
        'ivar consensus -p {params.prefix:q} -q {params.min_quality} -t {params.min_frequency} -m {params.min_depth} -c {params.min_ins_frequency} -n {params.char_under_min_depth} <{input.pileup:q} >{log:q} 2>&1 && '
        'mv "{params.prefix}.fa" {output.fasta:q} && mv "{params.prefix}.qual.txt" {output.quality:q}'


rule variant_calling:
    group: "sample"
    conda: "../envs/reads.yaml"
    shadow: "minimal"
    params:
        min_quality = 20,
        min_frequency = 0.05,
        min_depth = 30
    input:
        reference = "data/snpEff/data/{}/sequences.fa".format(config["REFERENCE"]),
        pileup = OUTPUT/"variants/pileup/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.pileup"
    output:
        tsv = OUTPUT/"variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.tsv"
    resources:
        runtime = lambda wc, attempt: 20 * attempt,
        mem_mb = lambda wc, attempt: 4000 * attempt
    retries: 2
    log: OUTPUT/"logs/variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/variant_calling.txt"
    shell: "ivar variants -p result -q {params.min_quality} -t {params.min_frequency} -m {params.min_depth} -r {input.reference:q} <{input.pileup:q} >{log:q} 2>&1 && mv result.tsv {output.tsv:q}"

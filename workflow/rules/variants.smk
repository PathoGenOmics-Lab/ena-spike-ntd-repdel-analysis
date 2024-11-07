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
    resources:
        runtime = lambda wc, attempt: 15 * attempt,
        mem_gb = lambda wc, attempt: 4 * attempt
    retries: 2
    log: "output/logs/variants/pileup/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}.txt"
    shell: "samtools mpileup -aa -x -A -d {params.max_depth} -B -Q {params.min_quality} -f {input.reference:q} {input.bam:q} >{output.pileup:q} 2>{log:q}"


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
        char_under_min_depth = "N",
        prefix = "{study}__{sample}__{platform}__{run}__{layout}__{strategy}"
    output:
        fasta = "output/variants/consensus/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.fasta",
        quality = "output/variants/consensus/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.qual.txt"
    log: "output/logs/variants/consensus/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}.txt"
    resources:
        runtime = lambda wc, attempt: 15 * attempt,
        mem_gb = lambda wc, attempt: 4 * attempt
    retries: 2
    shell:
        'ivar consensus -p {params.prefix:q} -q {params.min_quality} -t {params.min_frequency} -m {params.min_depth} -c {params.min_ins_frequency} -n {params.char_under_min_depth} <{input.pileup:q} >{log:q} 2>&1 && '
        'mv "{params.prefix}.fa" {output.fasta:q} && mv "{params.prefix}.qual.txt" {output.quality:q}'


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
    resources:
        runtime = lambda wc, attempt: 15 * attempt,
        mem_gb = lambda wc, attempt: 4 * attempt
    retries: 2
    log: "output/logs/variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/variant_calling.txt"
    shell: "ivar variants -p result -q {params.min_quality} -t {params.min_frequency} -m {params.min_depth} -g {input.gff} -r {input.reference:q} <{input.pileup:q} >{log:q} 2>&1 && mv result.tsv {output.tsv:q}"


rule ivar_tsv_to_vcf:
    conda: "../envs/pydata.yaml"
    input:
        tsv = "output/variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.tsv",
        reference = "output/reference/sequence.fasta"
    params:
        pass_only = False,
        allele_freq_threshold = 0,
        ignore_strand_bias = False,
        ignore_merge_codons = False,
        sample_name = "{study}_{sample}_{platform}_{run}_{layout}_{strategy}"
    output:
        vcf = "output/variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.vcf"
    log: "output/logs/variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/ivar_tsv_to_vcf.txt"
    script: "../scripts/ivar_tsv_to_vcf.py"


rule snpeff_annotate:
    conda: "../envs/annotation.yaml"
    input:
        datadir = "output/reference/snpeff/NC_045512.2",
        vcf = "output/variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.vcf"
    params:
        reference = "NC_045512.2"
    output:
        vcf = "output/variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.annotated.vcf"
    log: "output/logs/variants/snpeff_annotate/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}.txt"
    shell: "snpEff eff -dataDir {input.datadir:q} {params.reference} {input.vcf:q} >{output.vcf:q} 2>{log:q}"

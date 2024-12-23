rule ivar_tsv_to_vcf:
    group: "sample"
    conda: "../envs/pydata.yaml"
    input:
        tsv = OUTPUT/"variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.tsv",
        reference = OUTPUT/"reference/sequence.fasta"
    params:
        pass_only = False,
        allele_freq_threshold = 0,
        ignore_strand_bias = False,
        ignore_merge_codons = False,
        sample_name = "{study}_{sample}_{platform}_{run}_{layout}_{strategy}"
    output:
        vcf = temp(OUTPUT/"variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.vcf")
    resources:
        runtime = lambda wc, attempt: 5 * attempt,
        mem_mb = lambda wc, attempt: 100 * attempt
    retries: 2
    log: OUTPUT/"logs/variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/ivar_tsv_to_vcf.txt"
    script: "../scripts/ivar_tsv_to_vcf.py"


rule snpeff_annotate:
    group: "sample"
    conda: "../envs/annotation.yaml"
    shadow: "minimal"
    input:
        datadir = OUTPUT/"reference/snpeff/NC_045512.2",
        vcf = OUTPUT/"variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.vcf"
    params:
        reference = "NC_045512.2"
    output:
        vcf = temp(OUTPUT/"variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.annotated.vcf")
    resources:
        runtime = lambda wc, attempt: 5 * attempt,
        mem_mb = lambda wc, attempt: 200 * attempt
    retries: 2
    log: OUTPUT/"logs/variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/snpeff_annotate.txt"
    shell: "snpEff eff -dataDir {input.datadir:q} {params.reference} {input.vcf:q} >{output.vcf:q} 2>{log:q}"


rule snpsift_extract_variants:
    group: "sample"
    conda: "../envs/annotation.yaml"
    input:
        vcf = OUTPUT/"variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.annotated.vcf"
    params:
        min_depth = 40,
        hgvs_p_filter = build_snpsift_hgvs_p_filter,
        extract_columns = ["CHROM", "REF", "POS", "ALT", "DP", '"GEN[*].ALT_DP"', '"GEN[*].ALT_RV"', '"GEN[*].ALT_FREQ"', '"GEN[*].ALT_QUAL"', '"ANN[*].GENE"', '"ANN[*].HGVS_P"']
    output:
        tsv = temp(OUTPUT/"variants/snpsift_extract_variants/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/{haplotype}.tsv")
    resources:
        runtime = lambda wc, attempt: 1 * attempt,
        mem_mb = lambda wc, attempt: 200 * attempt
    retries: 2
    log:
        OUTPUT/"logs/variants/snpsift_extract_variants/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/{haplotype}_filter.txt",
        OUTPUT/"logs/variants/snpsift_extract_variants/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/{haplotype}_extractFields.txt"
    shell:
        'SnpSift filter "(DP >= {params.min_depth}){params.hgvs_p_filter}" {input.vcf:q} 2>{log[0]:q} | '
        'SnpSift extractFields -s "," - {params.extract_columns} >{output.tsv:q} 2>{log[1]:q}'

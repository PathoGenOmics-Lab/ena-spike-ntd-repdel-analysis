rule pileup:
    threads: 1
    group: "process"
    conda: "../envs/reads.yaml"
    input:
        reference = "output/reference/sequence.fasta",
        bam = "output/mapping/sorted_bam/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.sorted.bam"
    params:
        max_depth = 0,   # 0 means unrestricted
        min_quality = 0  # filtered later with iVar
    output:
        pileup = temp("output/variants/pileup/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.pileup")
    resources:
        runtime = lambda wc, attempt: 20 * attempt,
        mem_gb = lambda wc, attempt: 1 * attempt
    retries: 2
    log: "output/logs/variants/pileup/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}.txt"
    shell: "samtools mpileup -aa -x -A -d {params.max_depth} -B -Q {params.min_quality} -f {input.reference:q} {input.bam:q} >{output.pileup:q} 2>{log:q}"


rule coverage:
    threads: 1
    group: "process"
    shadow: "minimal"
    conda: "../envs/reads.yaml"
    input:
        bam = "output/mapping/sorted_bam/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.sorted.bam"
    params:
        chrom = "NC_045512.2",
        region_start = config["COVERAGE_FILTER"]["START"],
        region_end = config["COVERAGE_FILTER"]["END"]
    output:
        table = temp("output/variants/coverage/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.tsv")
    log: "output/logs/variants/coverage/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}.txt"
    shell:
        'printf "sample\n{wildcards.study}__{wildcards.sample}__{wildcards.platform}__{wildcards.run}__{wildcards.layout}__{wildcards.strategy}" >sample.txt && '
        "samtools coverage -d 0 -r {params.chrom}:{params.region_start}-{params.region_end} {input.bam:q} >coverage.txt 2>{log:q} && "
        "paste coverage.txt sample.txt >{output.table:q}"


rule coverage_merge:
    threads: 1
    shadow: "shallow"
    input: lambda w: build_search_targets(w, "output/variants/coverage/{}/{}/{}/{}/{}_{}_{}/sample.tsv")
    output: "output/variants/coverage.tsv"
    resources:
        runtime = "10m",
        mem_mb = 500
    shell: "head -n 1 {input[0]:q} >{output:q} && tail -n +2 -q {input:q} >>{output:q}"


rule filter_coverage:
    threads: 1
    input: "output/variants/coverage.tsv"
    params:
        # -1 == unfiltered
        min_threshold = {
            "numreads":  config["COVERAGE_FILTER"]["MIN_NUMREADS"],
            "covbases":  config["COVERAGE_FILTER"]["MIN_COVBASES"],
            "coverage":  config["COVERAGE_FILTER"]["MIN_COVERAGE"],
            "meandepth": config["COVERAGE_FILTER"]["MIN_MEANDEPTH"],
            "meanbaseq": config["COVERAGE_FILTER"]["MIN_MEANBASEQ"],
            "meanmapq":  config["COVERAGE_FILTER"]["MIN_MEANMAPQ"]
        }
    output: "output/variants/coverage.filtered.csv"
    run:
        import logging
        import csv
        logging.basicConfig(
            level=logging.INFO,
            format=config["PY_LOG_FMT"],
            filename=log[0]
        )
        n = 0
        ntotal = 0
        with open(input.table) as f, open(output.table, "w") as fw:
            reader = csv.DictReader(f, delimiter="\t")
            writer = csv.DictWriter(fw, fieldnames=reader.fieldnames)
            writer.writeheader()
            for row in reader:
                if all(float(row[colname]) <= value for colname, value in params.min_threshold.items() if value > -1):
                    writer.writerow(row)
                    n += 1
                ntotal += 1
        logging.info(f"Wrote {n} of {ntotal} records with all min_threshold={params.min_threshold}")


rule consensus:
    threads: 1
    group: "ivar"
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
        runtime = lambda wc, attempt: 30 * attempt,
        mem_gb = lambda wc, attempt: 4 * attempt
    retries: 2
    shell:
        'ivar consensus -p {params.prefix:q} -q {params.min_quality} -t {params.min_frequency} -m {params.min_depth} -c {params.min_ins_frequency} -n {params.char_under_min_depth} <{input.pileup:q} >{log:q} 2>&1 && '
        'mv "{params.prefix}.fa" {output.fasta:q} && mv "{params.prefix}.qual.txt" {output.quality:q}'


rule variant_calling:
    threads: 1
    group: "ivar"
    conda: "../envs/reads.yaml"
    shadow: "minimal"
    params:
        min_quality = 20,
        min_frequency = 0.05,
        min_depth = 30
    input:
        reference = "output/reference/sequence.fasta",
        pileup = "output/variants/pileup/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.pileup"
    output:
        tsv = "output/variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.tsv"
    resources:
        runtime = lambda wc, attempt: 20 * attempt,
        mem_gb = lambda wc, attempt: 4 * attempt
    retries: 2
    log: "output/logs/variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/variant_calling.txt"
    shell: "ivar variants -p result -q {params.min_quality} -t {params.min_frequency} -m {params.min_depth} -r {input.reference:q} <{input.pileup:q} >{log:q} 2>&1 && mv result.tsv {output.tsv:q}"


rule ivar_tsv_to_vcf:
    threads: 1
    group: "ivar"
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
        vcf = temp("output/variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.vcf")
    resources:
        runtime = lambda wc, attempt: 5 * attempt,
        mem_mb = lambda wc, attempt: 100 * attempt
    retries: 2
    log: "output/logs/variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/ivar_tsv_to_vcf.txt"
    script: "../scripts/ivar_tsv_to_vcf.py"


rule snpeff_annotate:
    threads: 1
    group: "snp"
    conda: "../envs/annotation.yaml"
    input:
        datadir = "output/reference/snpeff/NC_045512.2",
        vcf = "output/variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.vcf"
    params:
        reference = "NC_045512.2"
    output:
        vcf = temp("output/variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.annotated.vcf")
    resources:
        runtime = lambda wc, attempt: 5 * attempt,
        mem_mb = lambda wc, attempt: 200 * attempt
    retries: 2
    log: "output/logs/variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/snpeff_annotate.txt"
    shell: "snpEff eff -dataDir {input.datadir:q} {params.reference} {input.vcf:q} >{output.vcf:q} 2>{log:q}"


rule snpsift_extract_variants:
    threads: 1
    group: "snp"
    conda: "../envs/annotation.yaml"
    input:
        vcf = "output/variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.annotated.vcf"
    params:
        min_depth = 40,
        hgvs_p_filter = build_snpsift_hgvs_p_filter,
        extract_columns = ["CHROM", "REF", "POS", "ALT", "DP", '"GEN[*].ALT_DP"', '"GEN[*].ALT_RV"', '"GEN[*].ALT_FREQ"', '"GEN[*].ALT_QUAL"', '"ANN[*].GENE"', '"ANN[*].HGVS_P"']
    output:
        tsv = temp("output/variants/snpsift_extract_variants/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/{haplotype}.tsv")
    resources:
        runtime = lambda wc, attempt: 1 * attempt,
        mem_mb = lambda wc, attempt: 200 * attempt
    retries: 2
    log:
        "output/logs/variants/snpsift_extract_variants/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/{haplotype}_filter.txt",
        "output/logs/variants/snpsift_extract_variants/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/{haplotype}_extractFields.txt"
    shell:
        'SnpSift filter "(DP >= {params.min_depth}){params.hgvs_p_filter}" {input.vcf:q} 2>{log[0]:q} | '
        'SnpSift extractFields -s "," - {params.extract_columns} >{output.tsv:q} 2>{log[1]:q}'

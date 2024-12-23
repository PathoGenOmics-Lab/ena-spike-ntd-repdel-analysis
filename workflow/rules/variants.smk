rule pileup:
    group: "sample"
    conda: "../envs/reads.yaml"
    input:
        reference = OUTPUT/"reference/sequence.fasta",
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


rule coverage:
    group: "sample"
    shadow: "minimal"
    conda: "../envs/reads.yaml"
    input:
        bam = OUTPUT/"mapping/sorted_bam/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.sorted.bam"
    params:
        chrom = "NC_045512.2",
        region_start = config["COVERAGE_FILTER"]["START"],
        region_end = config["COVERAGE_FILTER"]["END"]
    output:
        table = temp(OUTPUT/"variants/coverage/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/coverage.tsv"),
        index = temp(OUTPUT/"mapping/sorted_bam/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.sorted.bam.bai")
    resources:
        runtime = lambda wc, attempt: 5 * attempt,
        mem_mb = lambda wc, attempt: 1000 * attempt
    retries: 2
    log: OUTPUT/"logs/variants/coverage/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}.txt"
    shell:
        'printf "sample\n{wildcards.study}__{wildcards.sample}__{wildcards.platform}__{wildcards.run}__{wildcards.layout}__{wildcards.strategy}" >sample.txt && '
        "samtools index --bai -o {output.index:q} {input.bam:q} && "
        "samtools coverage -d 0 -r {params.chrom}:{params.region_start}-{params.region_end} {input.bam:q} >coverage.txt 2>{log:q} && "
        "paste coverage.txt sample.txt >{output.table:q}"


rule filter_coverage:
    input: OUTPUT/"variants/coverage/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/coverage.tsv"
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
    output: OUTPUT/"variants/coverage/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/coverage.filtered.csv"
    log: OUTPUT/"logs/variants/coverage/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}.txt"
    run:
        import logging
        import csv
        logging.basicConfig(
            level=logging.INFO,
            format=config["PY_LOG_FMT"],
            filename=log[0]
        )
        passed = False
        with open(input[0]) as f, open(output[0], "w") as fw:
            reader = csv.DictReader(f, delimiter="\t")
            writer = csv.DictWriter(fw, fieldnames=reader.fieldnames)
            writer.writeheader()
            row = next(reader)
            if all(float(row[colname]) >= value for colname, value in params.min_threshold.items() if value > -1):
                writer.writerow(row)
                passed = True
        logging.info(f"Record pass={passed}")


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
        reference = OUTPUT/"reference/sequence.fasta",
        pileup = OUTPUT/"variants/pileup/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.pileup"
    output:
        tsv = OUTPUT/"variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/sample.tsv"
    resources:
        runtime = lambda wc, attempt: 20 * attempt,
        mem_mb = lambda wc, attempt: 4000 * attempt
    retries: 2
    log: OUTPUT/"logs/variants/variant_calling/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}/variant_calling.txt"
    shell: "ivar variants -p result -q {params.min_quality} -t {params.min_frequency} -m {params.min_depth} -r {input.reference:q} <{input.pileup:q} >{log:q} 2>&1 && mv result.tsv {output.tsv:q}"

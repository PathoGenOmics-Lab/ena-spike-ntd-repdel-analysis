rule fastqc:
    group: "fastqc"
    conda: "../envs/qc.yaml"
    shadow: "minimal"
    input: OUTPUT/"ena/downloads/fastq/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}"
    output: directory(OUTPUT/"preproc/fastqc/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}")
    log: OUTPUT/"logs/preproc/fastqc/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}.txt"
    shell: "mkdir out && fastqc --noextract -o out {input:q}/*.fastq.gz 2>{log} && mv out {output:q}"


rule multiqc:
    conda: "../envs/qc.yaml"
    input:
        lambda w: build_afterproc_targets_filtering(w, OUTPUT/f"preproc/fastp/{w.study}/{{}}/{{}}/{{}}/{{}}_{{}}_{{}}/report.json", ("sample_accession", "instrument_platform", "run_accession", "library_layout", "fastq_ftp", "library_strategy"), study_accession=w.study)
    output:
        directory(OUTPUT/"preproc/multiqc/{study}")
    log: OUTPUT/"logs/preproc/multiqc/{study}.txt"
    shell: "mkdir -p {output:q} && multiqc --outdir {output:q} --dirs -dd 5 {input:q} 2>{log:q}"


rule fastp_single:
    group: "fastp"
    conda: "../envs/qc.yaml"
    input:
        fastq = OUTPUT/"ena/downloads/fastq/{study}/{sample}/{platform}/{run}/{layout}_1_{strategy}/sample.fastq.gz"
    output:
        report = OUTPUT/"preproc/fastp/{study}/{sample}/{platform}/{run}/{layout}_1_{strategy}/report.html",
        json = OUTPUT/"preproc/fastp/{study}/{sample}/{platform}/{run}/{layout}_1_{strategy}/report.json",
        fastq = temp(OUTPUT/"preproc/fastq/{study}/{sample}/{platform}/{run}/{layout}_1_{strategy}/sample.fastp.fastq.gz")
    resources:
        runtime = lambda wc, attempt: 15 * attempt,
        mem_gb = lambda wc, attempt: 4 * attempt
    retries: 2
    log: OUTPUT/"logs/preproc/fastp_single/{study}/{sample}/{platform}/{run}/{layout}_1_{strategy}.txt"
    shell: "fastp -i {input.fastq:q} -o {output.fastq:q} -h {output.report:q} -j {output.json:q} 2>{log}"


rule fastp_paired:
    group: "fastp"
    conda: "../envs/qc.yaml"
    input:
        fastq_1 = OUTPUT/"ena/downloads/fastq/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/sample.R1.fastq.gz",
        fastq_2 = OUTPUT/"ena/downloads/fastq/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/sample.R1.fastq.gz"
    output:
        report = OUTPUT/"preproc/fastp/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/report.html",
        json = OUTPUT/"preproc/fastp/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/report.json",
        fastq_1 = temp(OUTPUT/"preproc/fastq/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/sample.fastp.R1.fastq.gz"),
        fastq_2 = temp(OUTPUT/"preproc/fastq/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/sample.fastp.R2.fastq.gz")
    resources:
        runtime = lambda wc, attempt: 15 * attempt,
        mem_gb = lambda wc, attempt: 4 * attempt
    retries: 2
    log: OUTPUT/"logs/preproc/fastp_paired/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}.txt"
    shell: "fastp -i {input.fastq_1:q} -I {input.fastq_2:q} -o {output.fastq_1:q} -O {output.fastq_2:q} -h {output.report:q} -j {output.json:q} 2>{log}"

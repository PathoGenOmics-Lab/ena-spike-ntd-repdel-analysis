rule fastqc:
    threads: 1
    group: "process"
    conda: "../envs/qc.yaml"
    shadow: "minimal"
    input: "output/ena/downloads/fastq/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}"
    output: directory("output/preproc/fastqc/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}")
    log: "output/logs/preproc/fastqc/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}.txt"
    shell: "mkdir out && fastqc --noextract -o out {input:q}/*.fastq.gz 2>{log} && mv out {output:q}"


rule multiqc:
    conda: "../envs/qc.yaml"
    input:
        lambda w: build_pangolin_targets_filtering(w, f"output/preproc/fastp/{w.study}/{{}}/{{}}/{{}}/{{}}_{{}}_{{}}/report.json", ("sample_accession", "instrument_platform", "run_accession", "library_layout", "fastq_ftp", "library_strategy"), study_accession=w.study)
    output:
        directory("output/preproc/multiqc/{study}")
    log: "output/logs/preproc/multiqc/{study}.txt"
    shell: "mkdir -p {output:q} && multiqc --outdir {output:q} --dirs -dd 5 {input:q} 2>{log:q}"


rule fastp_single:
    threads: 1
    group: "process"
    conda: "../envs/qc.yaml"
    input:
        fastq = "output/ena/downloads/fastq/{study}/{sample}/{platform}/{run}/{layout}_1_{strategy}/sample.fastq.gz"
    output:
        report = "output/preproc/fastp/{study}/{sample}/{platform}/{run}/{layout}_1_{strategy}/report.html",
        json = "output/preproc/fastp/{study}/{sample}/{platform}/{run}/{layout}_1_{strategy}/report.json",
        fastq = temp("output/preproc/fastq/{study}/{sample}/{platform}/{run}/{layout}_1_{strategy}/sample.fastp.fastq.gz")
    resources:
        runtime = lambda wc, attempt: 15 * attempt,
        mem_gb = lambda wc, attempt: 4 * attempt
    retries: 2
    log: "output/logs/preproc/fastp_single/{study}/{sample}/{platform}/{run}/{layout}_1_{strategy}.txt"
    shell: "fastp -i {input.fastq:q} -o {output.fastq:q} -h {output.report:q} -j {output.json:q} 2>{log}"


rule fastp_paired:
    threads: 1
    group: "process"
    conda: "../envs/qc.yaml"
    input:
        fastq_1 = "output/ena/downloads/fastq/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/sample.R1.fastq.gz",
        fastq_2 = "output/ena/downloads/fastq/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/sample.R1.fastq.gz"
    output:
        report = "output/preproc/fastp/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/report.html",
        json = "output/preproc/fastp/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/report.json",
        fastq_1 = temp("output/preproc/fastq/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/sample.fastp.R1.fastq.gz"),
        fastq_2 = temp("output/preproc/fastq/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/sample.fastp.R2.fastq.gz")
    resources:
        runtime = lambda wc, attempt: 15 * attempt,
        mem_gb = lambda wc, attempt: 4 * attempt
    retries: 2
    log: "output/logs/preproc/fastp_paired/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}.txt"
    shell: "fastp -i {input.fastq_1:q} -I {input.fastq_2:q} -o {output.fastq_1:q} -O {output.fastq_2:q} -h {output.report:q} -j {output.json:q} 2>{log}"

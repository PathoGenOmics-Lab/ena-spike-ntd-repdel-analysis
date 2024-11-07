rule fastqc:
    group: "preproc"
    conda: "../envs/qc.yaml"
    shadow: "minimal"
    input: "output/ena/downloads/fastq/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}"
    output: directory("output/preproc/fastqc/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}")
    log: "output/logs/preproc/fastqc/{study}/{sample}/{platform}/{run}/{layout}_{nfastq}_{strategy}.txt"
    shell: "mkdir out && fastqc --noextract -o out {input:q}/*.fastq.gz 2>{log} && mv out {output:q}"


rule multiqc:
    group: "preproc"
    conda: "../envs/qc.yaml"
    shadow: "minimal"
    input:
        lambda w: build_pangolin_targets(w, "output/preproc/fastp/{}/{}/{}/{}/{}_{}_{}/report.json")
    output:
        directory("output/preproc/multiqc/{study}")
    log: "output/logs/preproc/multiqc/{study}.txt"
    shell: "mkdir -p {output:q} && multiqc --outdir {output:q} --dirs -dd 5 {input:q} 2>{log:q}"


rule fastp_single:
    group: "preproc"
    conda: "../envs/qc.yaml"
    input:
        fastq = "output/ena/downloads/fastq/{study}/{sample}/{platform}/{run}/{layout}_1_{strategy}/sample.fastq.gz"
    output:
        report = "output/preproc/fastp/{study}/{sample}/{platform}/{run}/{layout}_1_{strategy}/report.html",
        json = "output/preproc/fastp/{study}/{sample}/{platform}/{run}/{layout}_1_{strategy}/report.json",
        fastq =temp("output/preproc/fastq/{study}/{sample}/{platform}/{run}/{layout}_1_{strategy}/sample.fastp.fastq.gz")
    log: "output/logs/preproc/fastp_single/{study}/{sample}/{platform}/{run}/{layout}_1_{strategy}.txt"
    shell: "fastp -i {input.fastq:q} -o {output.fastq:q} -h {output.report:q} -j {output.json:q} 2>{log}"


rule fastp_paired:
    group: "preproc"
    conda: "../envs/qc.yaml"
    input:
        fastq_1 = "output/ena/downloads/fastq/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/sample.R1.fastq.gz",
        fastq_2 = "output/ena/downloads/fastq/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/sample.R1.fastq.gz"
    output:
        report = "output/preproc/fastp/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/report.html",
        json = "output/preproc/fastp/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/report.json",
        fastq_1 = temp("output/preproc/fastq/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/sample.fastp.R1.fastq.gz"),
        fastq_2 = temp("output/preproc/fastq/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}/sample.fastp.R2.fastq.gz")
    log: "output/logs/preproc/fastp_paired/{study}/{sample}/{platform}/{run}/{layout}_2_{strategy}.txt"
    shell: "fastp -i {input.fastq_1:q} -I {input.fastq_2:q} -o {output.fastq_1:q} -O {output.fastq_2:q} -h {output.report:q} -j {output.json:q} 2>{log}"

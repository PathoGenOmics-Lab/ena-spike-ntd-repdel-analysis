rule fastqc:
    group: "preproc"
    conda: "../envs/qc.yaml"
    shadow: "minimal"
    input: "output/ena/downloads/fastq/{study}/{sample}/{platform}/{run}/{layout}_{strategy}"
    output: directory("output/preproc/fastqc/{study}/{sample}/{platform}/{run}/{layout}_{strategy}")
    shell: "mkdir out && fastqc --noextract -o out {input:q}/*.fastq.gz && mv out {output:q}"


rule multiqc:
    group: "preproc"
    conda: "../envs/qc.yaml"
    shadow: "minimal"
    input:
        lambda w: build_targets(w, "output/preproc/fastqc/{study}/{sample}/{platform}/{run}/{layout}_{strategy}"),
        lambda w: build_targets(w, "output/preproc/fastp/{study}/{sample}/{platform}/{run}/{layout}_{strategy}/report.html")
    output:
        directory("output/preproc/multiqc/{study}")
    shell: "mkdir out && multiqc --outdir out --dirs -dd 5 {input:q} && mv out {output:q}"


rule fastp_single:
    group: "preproc"
    conda: "../envs/qc.yaml"
    input:
        folder = "output/ena/downloads/fastq/{study}/{sample}/{platform}/{run}/SINGLE_{strategy}"
    output:
        report = "output/preproc/fastp/{study}/{sample}/{platform}/{run}/SINGLE_{strategy}/report.html",
        json = "output/preproc/fastp/{study}/{sample}/{platform}/{run}/SINGLE_{strategy}/report.json",
        fastq = "output/preproc/fastq/{study}/{sample}/{platform}/{run}/SINGLE_{strategy}/sample.fastp.fastq.gz"
    shell: """fastqs=( {input.folder:q}/*.fastq.gz ) && fastp -i "${{fastqs[0]}}" -o {output.fastq:q} -h {output.report:q} -j {output.json:q}"""


rule fastp_paired:
    group: "preproc"
    conda: "../envs/qc.yaml"
    input:
        folder = "output/ena/downloads/fastq/{study}/{sample}/{platform}/{run}/PAIRED_{strategy}"
    output:
        report = "output/preproc/fastp/{study}/{sample}/{platform}/{run}/PAIRED_{strategy}/report.html",
        json = "output/preproc/fastp/{study}/{sample}/{platform}/{run}/PAIRED_{strategy}/report.json",
        fastq_1 = "output/preproc/fastq/{study}/{sample}/{platform}/{run}/PAIRED_{strategy}/sample.fastp.fastq.R1.gz",
        fastq_2 = "output/preproc/fastq/{study}/{sample}/{platform}/{run}/PAIRED_{strategy}/sample.fastp.fastq.R2.gz"
    shell: """fastqs=( {input.folder:q}/*.fastq.gz ) && fastp -i "${{fastqs[0]}}" -I "${{fastqs[1]}}" -o {output.fastq_1:q} -O {output.fastq_2:q} -h {output.report:q} -j {output.json:q}"""

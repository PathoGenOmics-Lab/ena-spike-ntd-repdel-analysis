rule fastqc:
    group: "preproc"
    conda: "../envs/qc.yaml"
    shadow: "minimal"
    input: "output/ena/downloads/fastq/{study}/{sample}/{platform}/{run}/{layout}_{strategy}"
    output: directory("output/preproc/fastqc/{study}/{sample}/{platform}/{run}/{layout}_{strategy}")
    log: "output/logs/preproc/fastqc/{study}/{sample}/{platform}/{run}/{layout}_{strategy}.txt"
    shell: "mkdir out && fastqc --noextract -o out {input:q}/*.fastq.gz 2>{log} && mv out {output:q}"


rule multiqc:
    group: "preproc"
    conda: "../envs/qc.yaml"
    shadow: "minimal"
    input:
        lambda w: build_targets(w, "output/preproc/fastqc/{}/{}/{}/{}/{}_{}"),
        lambda w: build_targets(w, "output/preproc/fastp/{}/{}/{}/{}/{}_{}/report.html")
    output:
        directory("output/preproc/multiqc/{study}")
    log: "output/logs/preproc/multiqc/{study}.txt"
    shell: "mkdir out && multiqc --outdir out --dirs -dd 5 {input:q} 2>{log} && mv out {output:q}"


rule fastp_single:
    group: "preproc"
    conda: "../envs/qc.yaml"
    input:
        folder = "output/ena/downloads/fastq/{study}/{sample}/{platform}/{run}/SINGLE_{strategy}"
    output:
        report = "output/preproc/fastp/{study}/{sample}/{platform}/{run}/SINGLE_{strategy}/report.html",
        json = "output/preproc/fastp/{study}/{sample}/{platform}/{run}/SINGLE_{strategy}/report.json",
        fastq = "output/preproc/fastq/{study}/{sample}/{platform}/{run}/SINGLE_{strategy}/sample.fastp.fastq.gz"
    log: "output/logs/preproc/fastp_single/{study}/{sample}/{platform}/{run}/SINGLE_{strategy}.txt"
    shell:
        """
        fastqs=( {input.folder:q}/*.fastq.gz )
        if [ ${{#fastqs[@]}} -eq 1 ]; then
            fastp -i "${{fastqs[0]}}" -o {output.fastq:q} -h {output.report:q} -j {output.json:q} 2>{log}
        else
            echo There are ${{#fastqs[@]}} FASTQ files in the input directory (must be 1)
            exit 1
        fi
        """


rule fastp_paired:
    group: "preproc"
    conda: "../envs/qc.yaml"
    input:
        folder = "output/ena/downloads/fastq/{study}/{sample}/{platform}/{run}/PAIRED_{strategy}"
    output:
        report = "output/preproc/fastp/{study}/{sample}/{platform}/{run}/PAIRED_{strategy}/report.html",
        json = "output/preproc/fastp/{study}/{sample}/{platform}/{run}/PAIRED_{strategy}/report.json",
        fastq_1 = "output/preproc/fastq/{study}/{sample}/{platform}/{run}/PAIRED_{strategy}/sample.fastp.R1.fastq.gz",
        fastq_2 = "output/preproc/fastq/{study}/{sample}/{platform}/{run}/PAIRED_{strategy}/sample.fastp.R2.fastq.gz"
    log: "output/logs/preproc/fastp_paired/{study}/{sample}/{platform}/{run}/PAIRED_{strategy}.txt"
    shell:
        """
        fastqs=( {input.folder:q}/*.fastq.gz )
        if [ ${{#fastqs[@]}} -eq 1 ]; then
            fastp -i "${{fastqs[0]}}" -o {output.fastq:q} -h {output.report:q} -j {output.json:q} 2>{log}
        elif [ ${{#fastqs[@]}} -eq 2 ]; then
            fastq1=( {input.folder:q}/*.R1.fastq.gz ) && fastq2=( {input.folder:q}/*.R2.fastq.gz )
            fastp -i "${{fastq1[0]}}" -I "${{fastq2[0]}}" -o {output.fastq_1:q} -O {output.fastq_2:q} -h {output.report:q} -j {output.json:q} 2>{log}
        else
            echo There are ${{#fastqs[@]}} FASTQ files in the input directory (must be 1 or 2)
            exit 1
        fi
        """

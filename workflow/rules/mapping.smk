rule map_index:
    conda: "../envs/reads.yaml"
    input: "data/snpEff/data/{}/sequences.fa".format(config["REFERENCE"])
    output: OUTPUT/"reference/{preset}.mmi"
    resources:
        runtime = lambda wc, attempt: 15 * attempt,
        mem_mb = lambda wc, attempt: 8000 * attempt
    retries: 2
    log: OUTPUT/"logs/mapping/map_index/{preset}.txt"
    shell: "minimap2 -x {wildcards.preset:q} -d {output:q} {input:q} 2>{log}"


rule map_single_nanopore:
    group: "sample"
    threads: 2
    conda: "../envs/reads.yaml"
    input:
        reference = OUTPUT/"reference/map-ont.mmi",
        fastq = OUTPUT/"preproc/fastq/{study}/{sample}/OXFORD_NANOPORE/{run}/{layout}_1_{strategy}/sample.fastp.fastq.gz"
    output:
        bam = OUTPUT/"mapping/sorted_bam/{study}/{sample}/OXFORD_NANOPORE/{run}/{layout}_1_{strategy}/sample.sorted.bam"
    resources:
        runtime = lambda wc, attempt: 20 * attempt,
        mem_mb = lambda wc, attempt: 4000 * attempt
    retries: 2
    log:
        OUTPUT/"logs/mapping/map_single_nanopore/{study}/{sample}/OXFORD_NANOPORE/{run}/{layout}_1_{strategy}_minimap2.txt",
        OUTPUT/"logs/mapping/map_single_nanopore/{study}/{sample}/OXFORD_NANOPORE/{run}/{layout}_1_{strategy}_samtools.txt"
    shell: "minimap2 -t {threads} -ax map-ont {input.reference:q} {input.fastq:q} 2>{log[0]:q} | samtools sort -o {output.bam:q} 2>{log[1]:q}"


rule map_paired_illumina:
    group: "sample"
    threads: 2
    conda: "../envs/reads.yaml"
    input:
        reference = OUTPUT/"reference/sr.mmi",
        fastq_1 = OUTPUT/"preproc/fastq/{study}/{sample}/ILLUMINA/{run}/{layout}_2_{strategy}/sample.fastp.R1.fastq.gz",
        fastq_2 = OUTPUT/"preproc/fastq/{study}/{sample}/ILLUMINA/{run}/{layout}_2_{strategy}/sample.fastp.R2.fastq.gz"
    output:
        bam = OUTPUT/"mapping/sorted_bam/{study}/{sample}/ILLUMINA/{run}/{layout}_2_{strategy}/sample.sorted.bam"
    resources:
        runtime = lambda wc, attempt: 20 * attempt,
        mem_mb = lambda wc, attempt: 4000 * attempt
    retries: 2
    log:
        OUTPUT/"logs/mapping/map_paired_illumina/{study}/{sample}/ILLUMINA/{run}/{layout}_2_{strategy}_minimap2.txt",
        OUTPUT/"logs/mapping/map_paired_illumina/{study}/{sample}/ILLUMINA/{run}/{layout}_2_{strategy}_samtools.txt"
    shell: "minimap2 -t {threads} -ax sr {input.reference:q} {input.fastq_1:q} {input.fastq_2:q} 2>{log[0]:q} | samtools sort -o {output.bam:q} 2>{log[1]:q}"


rule map_single_illumina:
    group: "sample"
    threads: 2
    conda: "../envs/reads.yaml"
    input:
        reference = OUTPUT/"reference/sr.mmi",
        fastq = OUTPUT/"preproc/fastq/{study}/{sample}/ILLUMINA/{run}/{layout}_1_{strategy}/sample.fastp.fastq.gz"
    output:
        bam = OUTPUT/"mapping/sorted_bam/{study}/{sample}/ILLUMINA/{run}/{layout}_1_{strategy}/sample.sorted.bam"
    resources:
        runtime = lambda wc, attempt: 20 * attempt,
        mem_mb = lambda wc, attempt: 4000 * attempt
    retries: 2
    log:
        OUTPUT/"logs/mapping/map_single_illumina/{study}/{sample}/ILLUMINA/{run}/{layout}_1_{strategy}_minimap2.txt",
        OUTPUT/"logs/mapping/map_single_illumina/{study}/{sample}/ILLUMINA/{run}/{layout}_1_{strategy}_samtools.txt"
    shell: "minimap2 -t {threads} -ax sr {input.reference:q} {input.fastq:q} 2>{log[1]:q} | samtools sort -o {output.bam:q} 2>{log[1]:q}"


use rule map_single_illumina as map_single_ion_torrent with:
    group: "sample"
    conda: "../envs/reads.yaml"
    input:
        reference = OUTPUT/"reference/sr.mmi",
        fastq = OUTPUT/"preproc/fastq/{study}/{sample}/ION_TORRENT/{run}/{layout}_1_{strategy}/sample.fastp.fastq.gz"
    log:
        OUTPUT/"logs/mapping/map_single_ion_torrent/{study}/{sample}/ION_TORRENT/{run}/{layout}_1_{strategy}_minimap2.txt",
        OUTPUT/"logs/mapping/map_single_ion_torrent/{study}/{sample}/ION_TORRENT/{run}/{layout}_1_{strategy}_samtools.txt"
    output:
        bam = OUTPUT/"mapping/sorted_bam/{study}/{sample}/ION_TORRENT/{run}/{layout}_1_{strategy}/sample.sorted.bam"


use rule map_paired_illumina as map_paired_ion_torrent with:
    group: "sample"
    conda: "../envs/reads.yaml"
    input:
        reference = OUTPUT/"reference/sr.mmi",
        fastq_1 = OUTPUT/"preproc/fastq/{study}/{sample}/ION_TORRENT/{run}/{layout}_2_{strategy}/sample.fastp.R1.fastq.gz",
        fastq_2 = OUTPUT/"preproc/fastq/{study}/{sample}/ION_TORRENT/{run}/{layout}_2_{strategy}/sample.fastp.R2.fastq.gz"
    output:
        bam = OUTPUT/"mapping/sorted_bam/{study}/{sample}/ION_TORRENT/{run}/{layout}_2_{strategy}/sample.sorted.bam"
    log:
        OUTPUT/"logs/mapping/map_paired_ion_torrent/{study}/{sample}/ION_TORRENT/{run}/{layout}_2_{strategy}_minimap2.txt",
        OUTPUT/"logs/mapping/map_paired_ion_torrent/{study}/{sample}/ION_TORRENT/{run}/{layout}_2_{strategy}_samtools.txt"


rule map_pacbio_hifi:
    group: "sample"
    threads: 2
    conda: "../envs/reads.yaml"
    input:
        reference = OUTPUT/"reference/map-hifi.mmi",
        fastq = OUTPUT/"preproc/fastq/{study}/{sample}/PACBIO_SMRT/{run}/{layout}_1_{strategy}/sample.fastp.fastq.gz"
    output:
        bam = OUTPUT/"mapping/sorted_bam/{study}/{sample}/PACBIO_SMRT/{run}/{layout}_1_{strategy}/sample.sorted.bam"
    resources:
        runtime = lambda wc, attempt: 20 * attempt,
        mem_mb = lambda wc, attempt: 4000 * attempt
    retries: 2
    log:
        OUTPUT/"logs/mapping/map_pacbio_hifi/{study}/{sample}/PACBIO_SMRT/{run}/{layout}_1_{strategy}_minimap2.txt",
        OUTPUT/"logs/mapping/map_pacbio_hifi/{study}/{sample}/PACBIO_SMRT/{run}/{layout}_1_{strategy}_samtools.txt"
    shell: "minimap2 -t {threads} -ax map-hifi {input.reference:q} {input.fastq:q} 2>{log[0]:q} | samtools sort -o {output.bam:q} 2>{log[1]:q}"

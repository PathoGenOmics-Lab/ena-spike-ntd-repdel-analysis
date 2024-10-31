rule map_index:
    group: "mapping"
    conda: "../envs/reads.yaml"
    input: "output/reference/sequence.fasta"
    output: "output/reference/{preset}.mmi"
    resources:
        runtime = "15m"
    shell: "minimap2 -x {wildcards.preset:q} -d {output:q} {input:q}"


rule map_single_nanopore:
    threads: 4
    group: "mapping"
    conda: "../envs/reads.yaml"
    input:
        reference = "output/reference/map-ont.mmi",
        fastq = "output/preproc/fastq/{study}/{sample}/OXFORD_NANOPORE/{run}/SINGLE_{strategy}/sample.fastp.fastq.gz"
    output:
        bam = "output/mappings/sorted_bam/{study}/{sample}/OXFORD_NANOPORE/{run}/SINGLE_{strategy}/sample.sorted.bam"
    resources:
        runtime = "15m"
    shell: "minimap2 -t {threads} -ax map-ont {input.reference:q} {input.fastq:q} | samtools sort -o {output.bam:q}"


rule map_paired_illumina:
    threads: 4
    group: "mapping"
    conda: "../envs/reads.yaml"
    input:
        reference = "output/reference/sr.mmi",
        fastq_1 = "output/preproc/fastq/{study}/{sample}/ILLUMINA/{run}/PAIRED_{strategy}/sample.fastp.fastq.R1.gz",
        fastq_2 = "output/preproc/fastq/{study}/{sample}/ILLUMINA/{run}/PAIRED_{strategy}/sample.fastp.fastq.R2.gz"
    output:
        bam = "output/mappings/sorted_bam/{study}/{sample}/ILLUMINA/{run}/PAIRED_{strategy}/sample.sorted.bam"
    resources:
        runtime = "15m"
    shell: "minimap2 -t {threads} -ax sr {input.reference:q} {input.fastq_1:q} {input.fastq_2:q} | samtools sort -o {output.bam:q}"


rule map_single_illumina:
    threads: 4
    group: "mapping"
    conda: "../envs/reads.yaml"
    input:
        reference = "output/reference/sr.mmi",
        fastq = "output/preproc/fastq/{study}/{sample}/ILLUMINA/{run}/SINGLE_{strategy}/sample.fastp.fastq.gz"
    output:
        bam = "output/mappings/sorted_bam/{study}/{sample}/ILLUMINA/{run}/SINGLE_{strategy}/sample.sorted.bam"
    resources:
        runtime = "15m"
    shell: "minimap2 -t {threads} -ax sr {input.reference:q} {input.fastq:q} | samtools sort -o {output.bam:q}"


use rule map_single_illumina as map_ion_torrent with:
    threads: 4
    group: "mapping"
    conda: "../envs/reads.yaml"
    input:
        reference = "output/reference/sr.mmi",
        fastq = "output/preproc/fastq/{study}/{sample}/ION_TORRENT/{run}/{layout}_{strategy}/sample.fastp.fastq.gz"
    output:
        bam = "output/mappings/sorted_bam/{study}/{sample}/ION_TORRENT/{run}/{layout}_{strategy}/sample.sorted.bam"

rule map_pacbio_hifi:
    threads: 4
    group: "mapping"
    conda: "../envs/reads.yaml"
    input:
        reference = "output/reference/map-hifi.mmi",
        fastq = "output/preproc/fastq/{study}/{sample}/PACBIO_SMRT/{run}/{layout}_{strategy}/sample.fastp.fastq.gz"
    output:
        bam = "output/mappings/sorted_bam/{study}/{sample}/PACBIO_SMRT/{run}/{layout}_{strategy}/sample.sorted.bam"
    resources:
        runtime = "15m"
    shell: "minimap2 -t {threads} -ax map-hifi {input.reference:q} {input.fastq:q} | samtools sort -o {output.bam:q}"

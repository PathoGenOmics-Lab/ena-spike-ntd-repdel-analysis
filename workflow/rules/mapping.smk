rule map_index:
    group: "mapping"
    conda: "../envs/reads.yaml"
    input: "output/reference/sequence.fasta"
    output: "output/reference/{preset}.mmi"
    shell: "minimap2 -x {wildcards.preset:q} -d {output:q} {input:q}"


rule map_single_nanopore:
    threads: 8
    group: "mapping"
    conda: "../envs/reads.yaml"
    input:
        reference = "output/reference/map-ont.mmi",
        fastq = "output/preproc/fastq/{study}/{sample}/OXFORD_NANOPORE/{run}/SINGLE_{strategy}/sample.fastp.fastq.gz"
    output:
        bam = "output/mappings/sorted_bam/{study}/{sample}/OXFORD_NANOPORE/{run}/SINGLE_{strategy}/sample.sorted.bam"
    shell: "minimap2 -t {threads} -ax map-ont {input.reference:q} {input.fastq:q} | samtools sort -o {output.bam:q}"


rule map_paired_illumina:
    threads: 8
    group: "mapping"
    conda: "../envs/reads.yaml"
    input:
        reference = "output/reference/sr.mmi",
        fastq_1 = "output/preproc/fastq/{study}/{sample}/ILLUMINA/{run}/PAIRED_{strategy}/sample.fastp.fastq.R1.gz",
        fastq_2 = "output/preproc/fastq/{study}/{sample}/ILLUMINA/{run}/PAIRED_{strategy}/sample.fastp.fastq.R2.gz"
    output:
        bam = "output/mappings/sorted_bam/{study}/{sample}/ILLUMINA/{run}/PAIRED_{strategy}/sample.sorted.bam"
    shell: "minimap2 -t {threads} -ax sr {input.reference:q} {input.fastq_1:q} {input.fastq_2:q} | samtools sort -o {output.bam:q}"


rule map_single_illumina:
    threads: 8
    group: "mapping"
    conda: "../envs/reads.yaml"
    input:
        reference = "output/reference/sr.mmi",
        fastq = "output/preproc/fastq/{study}/{sample}/ILLUMINA/{run}/SINGLE_{strategy}/sample.fastp.fastq.gz"
    output:
        bam = "output/mappings/sorted_bam/{study}/{sample}/ILLUMINA/{run}/SINGLE_{strategy}/sample.sorted.bam"
    shell: "minimap2 -t {threads} -ax sr {input.reference:q} {input.fastq:q} | samtools sort -o {output.bam:q}"


use rule map_single_illumina as map_single_ion_torrent with:
    threads: 8
    group: "mapping"
    conda: "../envs/reads.yaml"
    input:
        reference = "output/reference/sr.mmi",
        fastq = "output/preproc/fastq/{study}/{sample}/ION_TORRENT/{run}/SINGLE_{strategy}/sample.fastp.fastq.gz"
    output:
        bam = "output/mappings/sorted_bam/{study}/{sample}/ION_TORRENT/{run}/SINGLE_{strategy}/sample.sorted.bam"

# TODO: PACBIO_SMRT

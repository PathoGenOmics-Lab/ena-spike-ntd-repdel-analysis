# ENA spike NTD repaired deletion analysis

[![PGO badge](https://img.shields.io/badge/PathoGenOmics-Lab-yellow.svg)](https://pathogenomics.github.io/)
[![DOI](https://img.shields.io/badge/Manuscript-under_review-387088.svg)]()
[![Release](https://img.shields.io/github/v/release/PathoGenOmics-Lab/ena-spike-ntd-repdel-analysis)](https://github.com/PathoGenOmics-Lab/ena-spike-ntd-repdel-analysis/releases)
[![Snakemake](https://img.shields.io/badge/Snakemake-8.25.3-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

A Snakemake workflow with associated scripts used for detecting spike NTD repaired deletions in SARS-CoV-2 Omicron BA.1 lineage reads. The workflow processes sequencing data retrieved from the ENA Portal API, performing quality filtering, read mapping, variant calling, and classification of the deletion repair genotype. This pipeline was developed as part of a larger study. The associated manuscript is currently under review.

Results generated with this pipeline are available via DOI: [10.20350/digitalCSIC/17032](https://doi.org/10.20350/digitalCSIC/17032). We ran [Snakemake v8.25.3](https://snakemake.readthedocs.io/en/v8.25.3/getting_started/installation.html) with [Python v3.12.7](https://www.python.org/downloads/release/python-3127/).

## Workflow summary

This pipeline fetches and processes SARS-CoV-2 "read run" ENA records with a sample collection date between 1 November 2021 and 1 August 2022, filtered for NCBI taxonomy code 2697049 (SARS-CoV-2) and *Homo sapiens* host, excluding sequencing platforms DNBseq, Element and capillary sequencing, and RNAseq, transcriptomic, metagenomic, and metatranscriptomic library strategies. Then, the following steps are run for each resulting record:

1. FASTQ retrieval via the ENA metadata FTP URLs.
2. Read preprocessing and quality filtering using [`fastp`](https://github.com/OpenGene/fastp) v0.23.4.
3. Read mapping with [`minimap2`](https://github.com/lh3/minimap2) v2.28 against a BA.1 reference genome (GenBank: OZ070629.1) using the recommended presets depending on the run sequencing platform.
4. Consensus genome generation with [`samtools`](https://github.com/samtools/samtools) v1.20 and [`iVar`](https://github.com/andersen-lab/ivar) v1.4.3.
5. Lineage assignment using [`pangolin`](https://github.com/cov-lineages/pangolin) v4.3.
6. Variant calling with [`iVar`](https://github.com/andersen-lab/ivar) v1.4.3, annotated with [`SnpEff`](https://github.com/pcingola/SnpEff) v5.2 and filtered with [`SnpSift`](https://github.com/pcingola/SnpSift) v5.2.
7. Classification in three "haplotypes", based on the presence or absence of S gene deletions ΔH69/V70 and ΔV143/Y145. Alleles are encoded as insertions in [HGVS nomenclature](https://hgvs-nomenclature.org), given the reference genome:
   - `Rep_69_70`: repair of S:ΔH69/V70 (`S:p.Val67_Ile68dup` detected, `S:p.Asp140_His141insValTyrTyr` absent).
   - `Rep_143_145`: repair of S:ΔV143/Y145 (`S:p.Asp140_His141insValTyrTyr` detected, `S:p.Val67_Ile68dup` absent).
   - `Rep_Both`: repair of both deletions (`S:p.Val67_Ile68dup` and `S:p.Asp140_His141insValTyrTyr` detected).
8. Data summarization and visualization using [`ape`](https://github.com/emmanuelparadis/ape) v5.8, [`Rsamtools`](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html) v2.18.0, [`tidyverse`](https://github.com/tidyverse/tidyverse) v2.0.0, and [`ggpubr`](https://github.com/kassambara/ggpubr) v0.6.0 in R v4.3.3.

## Usage

This repository contains a Snakemake workflow for processing sequencing data from FASTQ retrieval to classification and result summarization. The pipeline is conceptualized in two main sections: (1) an independent, linear processing pipeline for each record, and (2) summarization tasks that aggregate results and generate reports. Due to the large dataset size, a `LIGHT` configuration flag is available to execute only the first section of the DAG, reducing computational load.

### 1. Data retrieval and chunking

- [`00a_run_search.sh`](/00a_run_search.sh): queries the ENA Portal API to retrieve sequencing records.
- [`00b_generate_chunks.sh`](/00b_generate_chunks.sh): splits survey results into manageable chunks for processing via SLURM job arrays. The dataset is divided into 16 groups, each containing up to 5000 chunks, with each chunk holding 16 records. This approach addressed Snakemake limitations when handling large DAGs at the time of execution. Chunk settings were set considering [our HPC](https://garnatxadoc.uv.es) resource limits.

### 2. Haplotype classification

- [`01_run_haplotypes_array_chunked.sh <group number>`](/01_run_haplotypes_array_chunked.sh): runs the analysis for a specified chunk group. Each execution launches up to 5000 SLURM jobs through a job array. This step must be run for each group. For the final manuscript, only a subset of groups was analyzed. This step is executed with `LIGHT=True`. Most parameters can be tweaked via [the Snakemake `config.yaml` file](/config/config.yaml).

### 3. Summary and reporting

- [`02_run_complete.sh`](/02_run_complete.sh): Executes the full workflow to generate summary tables, visualizations, and reports. This step is executed with `LIGHT=False`. Given computational constraints, the final manuscript analysis used [`scripts/summarize_results.py`](/scripts/summarize_results.py) instead, which parses result files to produce a summary table with key measurements, and data visualizations were created manually using the same code integrated into the workflow.

## Citation

The manuscript is currently under review.

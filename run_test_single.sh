#!/usr/bin/env bash
#SBATCH --job-name repdel
#SBATCH --mem 64GB
#SBATCH --cpus-per-task 32
#SBATCH --ntasks 1
#SBATCH --qos short
#SBATCH --time 12:00:00

snakemake --workflow-profile profiles/default --resources mem_gb=64 --cores 32 --config SEARCH_LIMIT=100

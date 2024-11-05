#!/usr/bin/env bash
#SBATCH --job-name repdel
#SBATCH --mem 96GB
#SBATCH --cpus-per-task 32
#SBATCH --ntasks 1
#SBATCH --qos short
#SBATCH --time 1-00:00:00

snakemake --workflow-profile profiles/default --resources mem_gb=96 --cores 32

#!/usr/bin/env bash
#SBATCH --job-name 1trepdel
#SBATCH --mem 64GB
#SBATCH --cpus-per-task 32
#SBATCH --ntasks 1
#SBATCH --qos short
#SBATCH --time 12:00:00
#SBATCH --output slurm-%x-%J.out

snakemake --workflow-profile profiles/default --resources mem_gb=64 --cores 32 --config SEARCH_LIMIT=1000 --keep-going

#!/usr/bin/env bash
#SBATCH --job-name strepdel
#SBATCH --mem 8GB
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 1
#SBATCH --qos short
#SBATCH --time 12:00:00
#SBATCH --output slurm-%x-%J.out

snakemake --workflow-profile profiles/garnatxa --config SEARCH_LIMIT=1000 --keep-going

#!/usr/bin/env bash
#SBATCH --job-name 1repdel
#SBATCH --mem 96GB
#SBATCH --cpus-per-task 32
#SBATCH --ntasks 1
#SBATCH --qos short
#SBATCH --time 1-00:00:00
#SBATCH --output slurm-%x-%j-%J.out

snakemake --workflow-profile profiles/default --resources mem_gb=96 --cores 32 --scheduler greedy

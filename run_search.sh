#!/usr/bin/env bash
#SBATCH --job-name search
#SBATCH --mem 16GB
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 1
#SBATCH --qos interactive
#SBATCH --part interactive
#SBATCH --time 1-00:00:00
#SBATCH --output slurm-%x-%J.out

snakemake --workflow-profile profiles/garnatxa --until summarize_ena_search

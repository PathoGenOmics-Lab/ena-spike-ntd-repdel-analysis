#!/usr/bin/env bash
#SBATCH --job-name 02-srepdel
#SBATCH --mem 80GB
#SBATCH --cpus-per-task 32
#SBATCH --ntasks 1
#SBATCH --qos short
#SBATCH --time 06:00:00
#SBATCH --output slurm-%x-%j.out

echo ">>> LAUNCHING [03-final]"
srun snakemake \
    --keep-going \
    --workflow-profile profiles/default \
    --config SEARCH_TABLE=search.shuffled.filtered.tsv

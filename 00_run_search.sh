#!/usr/bin/env bash
#SBATCH --job-name srepdel
#SBATCH --mem 16GB
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 1
#SBATCH --qos short
#SBATCH --time 1-00:00:00
#SBATCH --output slurm-%x-%J.out

echo ">>> SEARCHING"
srun python scripts/search_ena.py search.raw.complete.tsv --start-date "2021-11-01" --end-date "2022-08-01"
srun python scripts/filter_search_ena.py search.raw.complete.tsv search.filtered.complete.tsv
srun python scripts/subsample.py search.filtered.complete.tsv search.tsv --subsample 100000

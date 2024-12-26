#!/usr/bin/env bash
#SBATCH --job-name 00-srepdel
#SBATCH --mem 16GB
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 1
#SBATCH --qos short
#SBATCH --time 1-00:00:00
#SBATCH --output slurm-%x-%J.out

set -e

echo ">>> SEARCHING"
if [ ! -f search.raw.complete.tsv ]; then
    srun python scripts/search_ena.py search.raw.complete.tsv --start-date "2021-11-01" --end-date "2022-08-01"
fi

echo ">>> SHUFFLING"
if [ ! -f search.shuffled.complete.tsv ]; then
    srun python scripts/shuffle_search.py search.raw.complete.tsv search.shuffled.complete.tsv
fi

echo ">>> FILTERING"
if [ ! -f search.shuffled.filtered.tsv ]; then
srun python scripts/filter_urls.py <search.shuffled.complete.tsv >search.shuffled.filtered.tsv
fi

echo ">>> BUILDING DATABASE"
if [ ! -f data/ena.sqlite ]; then
srun python scripts/to_sqlite.py search.shuffled.filtered.tsv search.sqlite
fi

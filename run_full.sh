#!/usr/bin/env bash
#SBATCH --job-name repdel
#SBATCH --mem 16GB
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 1
#SBATCH --qos medium
#SBATCH --time 7-00:00:00
#SBATCH --output slurm-%x-%J.out

set -e

if [ ! -f search.tsv ]; then
    echo ">>> SEARCHING"
    srun python scripts/search_ena.py search.raw.complete.tsv --start-date "2021-11-01" --end-date "2022-08-01"
    srun python scripts/filter_search_ena.py search.raw.complete.tsv search.filtered.complete.tsv
    srun python scripts/subsample.py search.filtered.complete.tsv search.tsv --subsample 100000
fi

echo ">>> RUNNING SECTION [1]"
srun snakemake batched.done --workflow-profile profiles/garnatxa --config UNTIL_FILTER=True
echo ">>> RUNNING SECTION [2]"
srun snakemake search.filtered.tsv --workflow-profile profiles/garnatxa --config UNTIL_FILTER=True
echo ">>> RUNNING SECTION [3]"
srun snakemake --workflow-profile profiles/garnatxa --config UNTIL_FILTER=False

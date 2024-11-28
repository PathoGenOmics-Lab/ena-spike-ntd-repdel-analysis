#!/usr/bin/env bash
#SBATCH --job-name t_repdel
#SBATCH --mem 16GB
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 1
#SBATCH --qos short
#SBATCH --time 1-00:00:00
#SBATCH --output slurm-%x-%J.out

set -e

if [ ! -f search.tsv ]; then
    echo ">>> SEARCHING"
    python search_ena.py search.tsv --start-date "2021-11-01" --end-date "2022-08-01"
    python subsample.py search.complete.tsv search.tsv --subsample 100000
fi

echo ">>> RUNNING SECTION [1]"
srun snakemake --workflow-profile profiles/garnatxa --config UNTIL_FILTER=True batched.done
echo ">>> RUNNING SECTION [2]"
srun snakemake --workflow-profile profiles/garnatxa --config UNTIL_FILTER=True search.filtered.tsv
echo ">>> RUNNING SECTION [3]"
srun snakemake --workflow-profile profiles/garnatxa --config UNTIL_FILTER=False

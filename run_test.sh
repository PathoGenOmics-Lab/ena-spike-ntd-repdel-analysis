#!/usr/bin/env bash
#SBATCH --job-name t_repdel
#SBATCH --mem 16GB
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 1
#SBATCH --qos short
#SBATCH --time 1-00:00:00
#SBATCH --output slurm-%x-%J.out

set -e

NBATCHES=2
LIMIT=1000

if [ ! -f search.tsv ]; then
    echo ">>> SEARCHING"
    python search_ena.py search.tsv \
        --start-date "2021-11-01" --end-date "2022-08-01" --limit $LIMIT
fi

for i in $(seq 1 $NBATCHES); do
    echo ">>> [1] START BATCH $i of $NBATCHES"
    srun snakemake --workflow-profile profiles/garnatxa --config UNTIL_FILTER=True --batch batcher=$i/$NBATCHES batched.done
done

for i in $(seq 1 $NBATCHES); do
    echo ">>> [2] START BATCH $i of $NBATCHES"
    srun snakemake --workflow-profile profiles/garnatxa --config UNTIL_FILTER=False --batch batcher=$i/$NBATCHES
done

#!/usr/bin/env bash
#SBATCH --job-name 02-srepdel
#SBATCH --mem 16GB
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 1
#SBATCH --qos medium
#SBATCH --time 7-00:00:00
#SBATCH --output slurm-%x-%J.out

set -e

NBATCHES=50

if [ -z "$1" ] && [ -z "$2" ]; then
    FROM=1
    TO=$NBATCHES
else
    FROM=$1
    TO=$2
fi

for i in $(seq $FROM $TO); do
    echo ">>> LAUNCHING BATCH $i of $NBATCHES"
    srun snakemake search.filtered.tsv --workflow-profile profiles/garnatxa --config UNTIL_FILTER=True --batch batcher=$i/$NBATCHES
done

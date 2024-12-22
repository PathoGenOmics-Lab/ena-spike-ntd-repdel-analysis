#!/usr/bin/env bash
#SBATCH --job-name 02a-srepdel
#SBATCH --mem 16GB
#SBATCH --cpus-per-task 32
#SBATCH --ntasks 1
#SBATCH --qos short
#SBATCH --time 06:00:00
#SBATCH --output slurm-%x-%j.out

NBATCHES=5000

echo ">>> LAUNCHING [02-$NBATCHES]"
srun snakemake \
    --keep-going \
    --workflow-profile profiles/default \
    --config UNTIL_FILTER=False SEARCH_TABLE=search.shuffled.filtered.tsv \
    --batch all=$NBATCHES/$NBATCHES

#!/usr/bin/env bash
#SBATCH --job-name 01b-srepdel
#SBATCH --mem 16GB
#SBATCH --cpus-per-task 32
#SBATCH --ntasks 1
#SBATCH --qos short
#SBATCH --time 06:00:00
#SBATCH --output slurm-%x-%A_%a.out

NBATCHES=...

echo ">>> LAUNCHING [$i]"
srun snakemake batched.done \
    --keep-going \
    --workflow-profile profiles/default \
    --config UNTIL_FILTER=True SEARCH_TABLE=search.shuffled.complete.tsv \
    --batch batcher=$NBATCHES/$NBATCHES

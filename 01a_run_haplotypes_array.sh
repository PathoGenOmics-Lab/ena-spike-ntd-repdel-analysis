#!/usr/bin/env bash
#SBATCH --job-name 01a-srepdel
#SBATCH --mem 80GB
#SBATCH --cpus-per-task 32
#SBATCH --ntasks 1
#SBATCH --qos short
#SBATCH --time 06:00:00
#SBATCH --output slurm-%x-%A_%a.out
#SBATCH --array 1-4999%200

NBATCHES=5000
i=$SLURM_ARRAY_TASK_ID

echo ">>> LAUNCHING [03-$i]"
srun snakemake \
    --keep-going \
    --workflow-profile profiles/default \
    --config SEARCH_TABLE=search.shuffled.filtered.tsv \
    --batch filter_haplotype=$i/$NBATCHES

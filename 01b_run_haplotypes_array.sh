#!/usr/bin/env bash
#SBATCH --job-name 01b-srepdel
#SBATCH --mem 80GB
#SBATCH --cpus-per-task 32
#SBATCH --ntasks 1
#SBATCH --qos short
#SBATCH --time 06:00:00
#SBATCH --output slurm-%x-%j.out

NBATCHES=5000
i=$NBATCHES

echo ">>> LAUNCHING [01-$i]"
srun snakemake batcher.done \
    --keep-going \
    --workflow-profile profiles/default \
    --batch batcher=$i/$NBATCHES

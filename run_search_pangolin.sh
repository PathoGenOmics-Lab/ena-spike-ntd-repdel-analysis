#!/usr/bin/env bash
#SBATCH --job-name sp_repdel
#SBATCH --mem 16GB
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 1
#SBATCH --qos interactive
#SBATCH --part interactive
#SBATCH --time 1-00:00:00
#SBATCH --output slurm-%x-%J.out

srun snakemake --workflow-profile profiles/garnatxa --until summarize_ena_search

NBATCHES=100
for i in $(seq 1 $NBATCHES); do
    srun snakemake --keep-going --workflow-profile profiles/garnatxa --until filter_search_ena_with_pangolin --batch $i/$N
done

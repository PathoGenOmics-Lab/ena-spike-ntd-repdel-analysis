#!/usr/bin/env bash
#SBATCH --job-name sp_repdel
#SBATCH --mem 16GB
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 2
#SBATCH --qos short
#SBATCH --time 1-00:00:00
#SBATCH --output slurm-%x-%J.out

NBATCHES=1000

if [ ! -f summarize_ena_search.done ]; then
    echo ">>> SEARCHING"
    srun snakemake --workflow-profile profiles/garnatxa --until summarize_ena_search
    touch summarize_ena_search.done
fi

for i in $(seq 1 $NBATCHES); do
    echo ">>> LAUNCHING BATCH $i"
    srun -n1 snakemake --keep-going --workflow-profile profiles/garnatxa --until summarize_ena_search_after_processing --batch coverage_merge=$i/$NBATCHES
done

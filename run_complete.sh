#!/usr/bin/env bash
#SBATCH --job-name srepdel
#SBATCH --mem 32GB
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 1
#SBATCH --qos short
#SBATCH --time 1-00:00:00
#SBATCH --output slurm-%x-%J.out

set -e

NBATCHES=1000

for i in $(seq 1 $NBATCHES); do
    echo ">>> START BATCH $i"
    srun snakemake --workflow-profile profiles/garnatxa --config UNTIL_FILTER=False --batch filter_haplotype=$i/$NBATCHES
end

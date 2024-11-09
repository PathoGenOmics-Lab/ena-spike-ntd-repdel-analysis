#!/usr/bin/env bash
#SBATCH --job-name repdel
#SBATCH --mem 32GB
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 1
#SBATCH --qos short
#SBATCH --time 1-00:00:00

set -e

NBATCHES=100

for i in {1..$NBATCHES}; do
    echo ">>> START BATCH $i"
    snakemake --workflow-profile profiles/garnatxa --scheduler greedy --batch all=$i/$NBATCHES
end

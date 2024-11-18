#!/usr/bin/env bash
#SBATCH --job-name sp_repdel
#SBATCH --mem 16GB
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 1
#SBATCH --qos short
#SBATCH --time 1-00:00:00
#SBATCH --output slurm-%x-%J.out

set -e

NBATCHES=2
LIMIT=1000

if [ ! -f search.test.tsv ]; then
    echo ">>> SEARCHING"
    python search_ena.py search.test.tsv \
        --start-date "2021-11-01" --end-date "2022-08-01" --limit $LIMIT
fi

for i in $(seq 1 $NBATCHES); do
    echo ">>> [1] START BATCH $i of $NBATCHES"
    srun snakemake --config OUTPUT=output_test SEARCH_TABLE=search.test.tsv FILTERED_TABLE=search.test.filtered.tsv \
        --workflow-profile profiles/garnatxa --config UNTIL_FILTER=True --batch consensus_merge=$i/$NBATCHES search.test.filtered.tsv
done

for i in $(seq 1 $NBATCHES); do
    echo ">>> [2] START BATCH $i of $NBATCHES"
    srun snakemake --config OUTPUT=output_test SEARCH_TABLE=search.test.tsv FILTERED_TABLE=search.test.filtered.tsv \
        --workflow-profile profiles/garnatxa --config UNTIL_FILTER=False --batch filter_haplotype=$i/$NBATCHES
end

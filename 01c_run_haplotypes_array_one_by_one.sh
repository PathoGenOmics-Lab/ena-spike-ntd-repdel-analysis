#!/usr/bin/env bash
#SBATCH --job-name 01c-srepdel
#SBATCH --mem 64GB
#SBATCH --cpus-per-task 32
#SBATCH --ntasks 1
#SBATCH --qos short
#SBATCH --time 1-00:00:00
#SBATCH --output slurm-%x-%A_%a.out

CHUNK_SIZE=50
i=$SLURM_ARRAY_TASK_ID
TABLE="search.shuffled.filtered.tsv"
TABLE_LINES=$(wc -l "$TABLE" | cut -f1 -d' ')
MAX_CHUNK=$(( $TABLE_LINES / $CHUNK_SIZE + 1 ))

head $TABLE >head.tsv

for chunk in $(seq 0 $MAX_CHUNK ); do
    for sample in $(scripts/iter_samples.py --chunk $chunk --size $CHUNK_SIZE --nrows $(( $TABLE_LINES - 1 )) <$TABLE); do
        srun snakemake \
            "output/repdel/filter_haplotype/$sample/"{Rep_69_70,Rep_143_145,Rep_Both}".inclpct_"{95,75}".exclpct_"{5,25}".csv" \
            --keep-going \
            --workflow-profile profiles/default \
            --config SEARCH_TABLE=head.tsv  # unused within workflow
    done
done

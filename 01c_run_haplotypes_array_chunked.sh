#!/usr/bin/env bash
#SBATCH --job-name 01c-srepdel
#SBATCH --mem-per-cpu 2GB
#SBATCH --cpus-per-task 2
#SBATCH --ntasks 16
#SBATCH --qos short
#SBATCH --time 1-00:00:00
#SBATCH --output slurm-%x-%j.out

CHUNK_SIZE=$SLURM_NTASKS
i=$SLURM_ARRAY_TASK_ID

TABLE="search.shuffled.filtered.tsv"
N_CHUNKS=$(scripts/iter_samples.py --size $CHUNK_SIZE --count <$TABLE)
MAX_CHUNK=$(( $N_CHUNKS - 1 ))

head $TABLE >head.tsv  # unused within workflow, saves memory

srun --ntasks 1 -c 2 snakemake --conda-create-envs-only --config SEARCH_TABLE=head.tsv
srun --ntasks 1 -c 2 snakemake "output/ena.sqlite" --config SEARCH_TABLE=$TABLE

for chunk in $(seq 0 $MAX_CHUNK); do
    echo "$(date) | >>> RUNNING FOR CHUNK $chunk OF $N_CHUNKS"
    for sample in $(scripts/iter_samples.py --chunk $chunk --size $CHUNK_SIZE <$TABLE); do
        echo "$(date) | >>> RUNNING FOR SAMPLE $sample"
        srun --ntasks 1 -c $SLURM_CPUS_PER_TASK --mem-per-cpu $SLURM_MEM_PER_CPU --output slurm-%x-%j_%s.out \
            snakemake \
            "output/repdel/filter_haplotype/$sample/"{Rep_69_70,Rep_143_145,Rep_Both}".inclpct_"{95,75}".exclpct_"{5,25}".csv" \
            --nolock \
            --keep-going \
            --workflow-profile profiles/default --cores $SLURM_CPUS_PER_TASK \
            --config SEARCH_TABLE=head.tsv &
    done
    wait
done

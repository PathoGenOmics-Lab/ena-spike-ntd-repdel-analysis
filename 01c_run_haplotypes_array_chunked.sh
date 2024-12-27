#!/usr/bin/env bash
#SBATCH --job-name 01c-srepdel
#SBATCH --mem-per-cpu 2GB
#SBATCH --cpus-per-task 16
#SBATCH --ntasks 1
#SBATCH --qos short
#SBATCH --time 01:00:00
#SBATCH --output logs/slurm-%x-%A_%a.out
#SBATCH --array 0-78141%200

TABLE="search.shuffled.filtered.tsv"
CHUNK="chunks/chunk_$SLURM_ARRAY_TASK_ID.txt"

echo "$(date) | >>> RUNNING FOR CHUNK $SLURM_ARRAY_TASK_ID ($CHUNK) OF $SLURM_ARRAY_TASK_MAX ON $SLURM_CPUS_PER_TASK CPUs"

samples=$(paste -s -d, $CHUNK)
paths="output/repdel/filter_haplotype/{$samples}/{Rep_69_70,Rep_143_145,Rep_Both}.inclpct_{95,75}.exclpct_{5,25}.csv"

srun --nodes 1 --ntasks 1 -c $SLURM_CPUS_PER_TASK --mem-per-cpu $SLURM_MEM_PER_CPU \
    snakemake \
    $(eval echo $paths) \
    --nolock \
    --keep-going \
    --workflow-profile profiles/default --cores $SLURM_CPUS_PER_TASK \
    --config SEARCH_TABLE=$TABLE LIGHT=True

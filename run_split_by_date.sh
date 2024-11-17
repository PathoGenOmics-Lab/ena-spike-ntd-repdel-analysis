#!/usr/bin/env bash
#SBATCH --job-name split-repdel
#SBATCH --mem 16GB
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 1
#SBATCH --qos medium
#SBATCH --time 3-00:00:00
#SBATCH --output slurm-%x-%J.out

N=1000  # approx

echo ">>> RUNNING SEARCH"
srun snakemake --workflow-profile profiles/garnatxa --until summarize_ena_search --config OUTPUT=output_search

echo ">>> SPLITTING DATES"
srun python split_dates.py -n $N output_search/ena/search.filtered.tsv dates.tsv

while read -r start end n; do
    echo ">>> RUNNING $n samples from $start to $end"
    srun snakemake --workflow-profile profiles/garnatxa --config OUTPUT=output_${n}_${start}_${end} START_DATE=$start END_DATE=$end
done < dates.tsv

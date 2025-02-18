#!/usr/bin/env bash
#SBATCH --job-name 00b-srepdel
#SBATCH --mem 8GB
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 1
#SBATCH --qos short
#SBATCH --time 1-00:00:00
#SBATCH --output slurm-%x-%J.out

N_PER_GROUP=$(( 5000 * 16 ))  # SLURM MaxArraySize * chunk size

header=$(mktemp)
records=$(mktemp)
groups_dir=$(mktemp -d)

echo "$(date) | Extracting records"
srun tail -n +2 search.shuffled.filtered.tsv >$records
srun head -n 1 search.shuffled.filtered.tsv >$header
echo "$(date) | Splitting table"
srun split -l $N_PER_GROUP -d $records $groups_dir/group_
echo "$(date) | Removing records file"
srun rm $records

for group in $groups_dir/group_*; do
    name=${group##*/}
    echo "$(date) | Chunking $name"
    srun cat $header $group | scripts/chunk_samples.py --size 16 chunks/$name
done

echo "$(date) | Removing groups directory"
srun rm $groups_dir/*
rmdir $groups_dir

#!/usr/bin/env bash
#SBATCH --job-name 00-srepdel
#SBATCH --mem 8GB
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 1
#SBATCH --qos short
#SBATCH --time 1-00:00:00
#SBATCH --output slurm-%x-%J.out

srun scripts/chunk_samples.py --size 16 chunks < search.shuffled.filtered.tsv

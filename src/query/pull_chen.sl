#!/bin/bash
#
#SBATCH --job-name=pull_chen
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=express
#SBATCH --time=02:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=1-20

module load R/4.2.2-foss-2022b

# Pull out GWAS variants from each METSIM summary stats

# FILENAME=$(ls /data/scratch/rtakei/projects/2023/metqtl_project/data/chen_metqtl/GCST*.tsv | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
FILENAME=$(cat /data/scratch/rtakei/projects/2023/metqtl_project/failed_process.txt | sed 's/.gz//g' | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
COMPOUND=$(echo ${FILENAME%%_build*} | sed 's/.*\///g')

mkdir -p results/metqtl/chen/${COMPOUND}/

Rscript /data/scratch/rtakei/projects/2023/metqtl_project/src/query/pull_chen.R ${FILENAME}


#!/bin/bash
#
#SBATCH --job-name=schlosser_plasma_query
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=express
#SBATCH --time=00:15:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=1-1296

module load R

# Pull out GWAS variants from each plasma metQTL summary stats

FILENAME=$(ls /data/scratch/rtakei/projects/2023/metqtl_project/data/schlosser_metqtl/plasma/GCST*tsv.gz | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
COMPOUND=$(echo ${FILENAME%%_build*} | sed 's/.*\///g')

mkdir -p results/metqtl/schlosser/plasma/${COMPOUND}/

Rscript /data/scratch/rtakei/projects/2023/metqtl_project/src/query/pull_plasma_metqtl.R ${FILENAME}


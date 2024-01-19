#!/bin/bash
#
#SBATCH --job-name=metsim_coloc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=express
#SBATCH --time=00:15:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=1-1391

module load R

# Pull out GWAS variants from each METSIM summary stats

FILENAME=$(ls /data/scratch/rtakei/projects/2023/metqtl_project/data/metsim_data/C*.tsv | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
COMPOUND=$(echo ${FILENAME%%.tsv} | sed 's/.*\///g')

mkdir -p results/coloc_res/metsim/{full,male,female}/

Rscript /data/scratch/rtakei/projects/2023/metqtl_project/src/coloc/run_coloc.metsim.R ${COMPOUND} full
Rscript /data/scratch/rtakei/projects/2023/metqtl_project/src/coloc/run_coloc.metsim.R ${COMPOUND} male
Rscript /data/scratch/rtakei/projects/2023/metqtl_project/src/coloc/run_coloc.metsim.R ${COMPOUND} female


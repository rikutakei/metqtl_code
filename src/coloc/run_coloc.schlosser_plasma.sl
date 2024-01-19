#!/bin/bash
#
#SBATCH --job-name=schlosser_plasma_coloc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=express
#SBATCH --time=00:30:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=1-1296

module load R

# Run coloc for a specific plasma metQTL

FILENAME=$(ls /data/scratch/rtakei/projects/2023/metqtl_project/data/schlosser_metqtl/plasma/GCST*tsv.gz | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
COMPOUND=$(echo ${FILENAME%%_build*} | sed 's/.*\///g')

mkdir -p results/coloc_res/schlosser/plasma/{full,female,male}/

Rscript /data/scratch/rtakei/projects/2023/metqtl_project/src/coloc/run_coloc.schlosser.R ${COMPOUND} full plasma
Rscript /data/scratch/rtakei/projects/2023/metqtl_project/src/coloc/run_coloc.schlosser.R ${COMPOUND} male plasma
Rscript /data/scratch/rtakei/projects/2023/metqtl_project/src/coloc/run_coloc.schlosser.R ${COMPOUND} female plasma


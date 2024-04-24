#!/bin/bash
#
#SBATCH --job-name=schlosser_urine_coloc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=express
#SBATCH --time=01:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
##SBATCH --array=1-1401
#SBATCH --array=57,72,73,74,75,78,79,80,81

module load R/4.2.0-foss-2022b

# Run coloc for a specific urine metQTL

FILENAME=$(ls /data/scratch/rtakei/projects/2023/metqtl_project/data/schlosser_metqtl/urine/GCST*tsv.gz | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
COMPOUND=$(echo ${FILENAME%%_build*} | sed 's/.*\///g')

mkdir -p results/coloc_res/schlosser/urine/{full,female,male}/

Rscript /data/scratch/rtakei/projects/2023/metqtl_project/src/coloc/run_coloc.schlosser.R ${COMPOUND} full urine
Rscript /data/scratch/rtakei/projects/2023/metqtl_project/src/coloc/run_coloc.schlosser.R ${COMPOUND} male urine
Rscript /data/scratch/rtakei/projects/2023/metqtl_project/src/coloc/run_coloc.schlosser.R ${COMPOUND} female urine


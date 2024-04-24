#!/bin/bash
#
#SBATCH --job-name=process_chen_raw
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=express
#SBATCH --time=00:10:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=1-1400

# Process the downloaded Chen summary stats

FILENAME=$(ls /data/scratch/rtakei/projects/2023/metqtl_project/data/chen_metqtl/*.tsv.gz | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)

bash /data/scratch/rtakei/projects/2023/metqtl_project/src/chen/process_raw.sh ${FILENAME}

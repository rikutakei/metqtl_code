#! /bin/bash
#SBATCH --job-name=schlosser_dl_plasma
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=express
#SBATCH --time=00:15:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=1-1296%20

# Download Schloser PLasma metQTL data

URL=$(cat /data/scratch/rtakei/projects/2023/metqtl_project/data/schlosser_metqtl/plasma_urls.txt | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)

wget -c -P data/schlosser_metqtl/plasma/ ${URL}


#! /bin/bash
#SBATCH --job-name=chen_dl_plasma
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=express
#SBATCH --time=00:20:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=1-1400%20

# Download Chen et al PLasma metQTL data

URL=$(cat /data/scratch/rtakei/projects/2023/metqtl_project/data/chen_metqtl/chen_urls.txt | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)

wget -c -P /data/scratch/rtakei/projects/2023/metqtl_project/data/chen_metqtl/ ${URL}


#!/bin/bash
#
#SBATCH --job-name=prep_metsim_var_gatk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=express
#SBATCH --time=00:10:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=1-1391

# Pull out relevant columns and foramt it for GATK build conversion

FILENAME=$(ls /data/scratch/rtakei/projects/2023/metqtl_project/data/metsim_data/C*.tsv | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)

cut -f1-4 ${FILENAME} | awk '{print $1, $2, $1"_"$2, $3, $4, ".", ".", "."}' | tail -n+2 | tr ' ' '\t' > ${FILENAME%%.tsv}.pre_vcf

cat /data/scratch/rtakei/projects/2023/metqtl_project/data/liftover/vcf_header.txt ${FILENAME%%.tsv}.pre_vcf > ${FILENAME%%.tsv}.vcf && rm ${FILENAME%%.tsv}.pre_vcf

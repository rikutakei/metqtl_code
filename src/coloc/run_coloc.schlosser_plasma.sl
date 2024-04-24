#!/bin/bash
#
#SBATCH --job-name=schlosser_plasma_coloc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=express
##SBATCH --time=00:30:00
#SBATCH --time=01:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
##SBATCH --array=1-1296
#SBATCH --array=7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,560,561,562,563,564,565,567,568,569,570,571,572,573,574,575,576,577,578,579,1112,1113,1114,1115,1116,1117,1118,1119,1120,1121,1122,1123,1124,1125,1126,1132,1133,1134,1135,1136,1137

module load R/4.2.0-foss-2022b

# Run coloc for a specific plasma metQTL

FILENAME=$(ls /data/scratch/rtakei/projects/2023/metqtl_project/data/schlosser_metqtl/plasma/GCST*tsv.gz | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
COMPOUND=$(echo ${FILENAME%%_build*} | sed 's/.*\///g')

mkdir -p results/coloc_res/schlosser/plasma/{full,female,male}/

Rscript /data/scratch/rtakei/projects/2023/metqtl_project/src/coloc/run_coloc.schlosser.R ${COMPOUND} full plasma
Rscript /data/scratch/rtakei/projects/2023/metqtl_project/src/coloc/run_coloc.schlosser.R ${COMPOUND} male plasma
Rscript /data/scratch/rtakei/projects/2023/metqtl_project/src/coloc/run_coloc.schlosser.R ${COMPOUND} female plasma


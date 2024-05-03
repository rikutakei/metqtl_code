#! /bin/bash
#SBATCH --job-name=clump_summary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=16G
#SBATCH --partition=express
#SBATCH --time=02:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=1-33

module load PLINK/1.90-foss-2016a

INPUT=$(ls /scratch/rtakei/projects/2023/metqtl_project/data/mr_dat/*clump.txt | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)
OUT=${INPUT%%.txt}_result
KGP=/data/project/merrimanlab/reference_files/1kgp_data/per_ancestry/EUR/plink/

# Pull out independent and significant variant (i.e. lead SNPs) from summary
# stats
plink --bfile ${KGP}/1KGP_EUR_wgs.all_samples.rsid.biallelic --clump ${INPUT} --clump-p1 5e-8 --clump-r2 0.01 --clump-kb 500 --out ${OUT}

# After clumping, calculate the LD between all of the lead variants to make
# sure they are independent
awk 'NR > 1 {print $3}' ${OUT}.clumped | grep -v '^$' > ${OUT}.lead.txt
plink --bfile ${KGP}/1KGP_EUR_wgs.all_samples.rsid.biallelic --extract ${OUT}.lead.txt --r2 inter-chr --ld-window-r2 0.01 --out ${OUT}


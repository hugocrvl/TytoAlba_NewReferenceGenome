#!/bin/bash

#SBATCH --job-name=1_Red.slurm
#SBATCH --output=1.1_Red_%A.out

#SBATCH --mail-type=ALL
#SBATCH --mail-user=hugo.corval@unil.ch

#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=01:00:00

WORKDIR=/work/FAC/FBM/DEE/jgoudet/barn_owl/hcorval/Chapter_NewReferenceGenome/9_Repeats
cd ${WORKDIR}

echo ""
echo "[$(date)] Start of the job"
echo ""

echo "Repeats detection along Fa haplotype"

mkdir -p 1.2.1_Red_Masks_Output_Fa
mkdir -p 1.3.1_Red_Repeats_Output_Fa

/work/FAC/FBM/DEE/jgoudet/barn_owl/tcumer/TOOLS/RED/bin/Red \
    -gnm /work/FAC/FBM/DEE/jgoudet/barn_owl/hcorval/Chapter_NewReferenceGenome/9_RepeatsRed/0.1_FaReference \
    -msk 1.2.1_Red_Masks_Output_Fa -rpt 1.3.1_Red_Repeats_Output_Fa

echo "Repeats detection along Mo haplotype"

mkdir -p 1.2.2_Red_Masks_Output_Mo
mkdir -p 1.3.2_Red_Repeats_Output_Mo

/work/FAC/FBM/DEE/jgoudet/barn_owl/tcumer/TOOLS/RED/bin/Red \
    -gnm /work/FAC/FBM/DEE/jgoudet/barn_owl/hcorval/Chapter_NewReferenceGenome/9_RepeatsRed/0.2_MoReference \
    -msk 1.2.2_Red_Masks_Output_Mo -rpt 1.3.2_Red_Repeats_Output_Mo

echo ""
echo "[$(date)] End of the job"
echo ""
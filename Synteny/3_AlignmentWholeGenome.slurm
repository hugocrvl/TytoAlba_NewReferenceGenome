#!/bin/bash

#SBATCH --job-name=3_AlignmentWholeGenome.slurm
#SBATCH --output=3.1_AlignmentWholeGenome_%A.out

#SBATCH --mail-type=ALL
#SBATCH --mail-user=hugo.corval@unil.ch

#SBATCH --cpus-per-task=5
#SBATCH --mem=20G
#SBATCH --time=01:00:00

echo ""
echo "[$(date)] Start of the job"
echo ""

module load minimap2

WORKDIR="/work/FAC/FBM/DEE/jgoudet/barn_owl/hcorval/Chapter_NewReferenceGenome/7_Synteny"
cd ${WORKDIR}

TARGET="/work/FAC/FBM/DEE/jgoudet/barn_owl/hcorval/Chapter_NewReferenceGenome/0_FinalFiles/Final_FatherAssembly.fa"
QUERY="/work/FAC/FBM/DEE/jgoudet/barn_owl/hcorval/Chapter_NewReferenceGenome/0_FinalFiles/Final_MotherAssembly.fa"
OUTPUT="./3.2_AlignmentWholeGenome.paf"

minimap2 -x asm10 -c --eqx --secondary=no ${TARGET} ${QUERY} > ${OUTPUT}

echo ""
echo "[$(date)] End of the job"
echo ""
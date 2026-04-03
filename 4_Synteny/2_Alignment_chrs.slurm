#!/bin/bash

#SBATCH --job-name=2_Alignment_chrs.slurm
#SBATCH --output=2.1_Alignment_chrs_%A.out

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

mkdir -p 2.2_Alignment_chrs

module load minimap2
for Fa_chr in $(cat 1.1_Chromosomes_Father.txt); do
    chr=$(echo "$Fa_chr" | sed 's/Fa_//')
    echo "${chr}"
    Fa_fasta=$(realpath 1.4_SplitAssembliesPerLGs_Father/Fa_${chr}.fasta)
    Mo_fasta=$(realpath 1.5_SplitAssembliesPerLGs_Mother/Mo_${chr}.fasta)
    output=2.2_Aligned_chrs/Aligned_${chr}.paf
    minimap2 -x asm10 -c --eqx --secondary=no ${Fa_fasta} ${Mo_fasta} > ${output}
done

echo ""
echo "[$(date)] End of the job"
echo ""
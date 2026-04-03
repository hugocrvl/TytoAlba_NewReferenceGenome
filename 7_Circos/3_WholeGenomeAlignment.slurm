#!/bin/bash

#SBATCH --job-name=3_WholeGenomeAlignment.slurm
#SBATCH --output=3.1_WholeGenomeAlignment_%A.out

#SBATCH --mail-type=ALL
#SBATCH --mail-user=hugo.corval@unil.ch

#SBATCH --cpus-per-task=5
#SBATCH --mem=20G
#SBATCH --time=00:01:00

if [ $# -ne 3 ]; then
  echo "Usage: $0 TARGET.fasta QUERY.fasta OUTPUT.paf" >&2
  exit 1
fi

TARGET=${1}
QUERY=${2}
OUTPUT=${3}

echo ""
echo "[$(date)] Start of the job"
echo ""

module load minimap2

WORKDIR="/work/FAC/FBM/DEE/jgoudet/barn_owl/hcorval/Chapter_NewReferenceGenome/11_Circos"
cd ${WORKDIR}

minimap2 -x asm5 -c --eqx -t ${SLURM_CPUS_PER_TASK} --secondary=no --cs=long ${TARGET} ${QUERY} > ${OUTPUT}

echo ""
echo "[$(date)] End of the job"
echo ""

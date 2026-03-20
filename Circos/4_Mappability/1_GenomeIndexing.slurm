#!/bin/bash

#SBATCH --job-name=1_GenomeIndexing.slurm
#SBATCH --output=4_Mappability/1.1_GenomeIndexing_%A.out

#SBATCH --mail-type=ALL
#SBATCH --mail-user=hugo.corval@unil.ch

#SBATCH --cpus-per-task=3
#SBATCH --mem=10G
#SBATCH --time=00:30:00

if [ $# -ne 2 ]; then
  echo "Usage: $0 FASTA OUTPUT_FILE" >&2
  exit 1
fi

FASTA=${1}
OUTPUT_DIR=${2}

echo ""
echo "[$(date)] Start of the job"
echo ""

module load miniforge3
conda_init
conda activate genmap

genmap index -F ${FASTA} -I ${OUTPUT_DIR}

echo ""
echo "[$(date)] End of the job"
echo ""

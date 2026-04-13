#!/bin/bash

#SBATCH --job-name=2_GenomeMappability.slurm
#SBATCH --output=4_Mappability/2.1_GenomeMappability_%A.out

#SBATCH --mail-type=ALL
#SBATCH --mail-user=hugo.corval@unil.ch

#SBATCH --cpus-per-task=5
#SBATCH --mem=20G
#SBATCH --time=48:00:00

if [ $# -ne 2 ]; then
  echo "Usage: $0 FASTA OUTPUT_FILE" >&2
  exit 1
fi

INDEX_DIR=${1}
OUTPUT_DIR=${2}

echo ""
echo "[$(date)] Start of the job"
echo ""

module load miniforge3
conda_init
conda activate genmap

mkdir -p ${OUTPUT_DIR}

genmap map -K 150 -E 2 -I ${INDEX_DIR} -O ${OUTPUT_DIR} -t -w -bg

echo ""
echo "[$(date)] End of the job"
echo ""
#!/bin/bash

#SBATCH --job-name=2.1_RunRepeatModeler.slurm
#SBATCH --output=2.1_RunRepeatModeler_%A.out

#SBATCH --mail-type=ALL
#SBATCH --mail-user=hugo.corval@unil.ch

#SBATCH --cpus-per-task=30
#SBATCH --mem=20G
#SBATCH --time=30:00:00

echo ""
echo "Start of the job [$(date)]"
echo ""

if [ $# -ne 2 ]; then
  echo "Usage: $0 FASTA DATABASE" >&2
  exit 1
fi

REF=${1}
DATABASE_DIR=${2}

cd ${DATABASE_DIR}

module load singularityce

export SINGULARITY_BINDPATH="/users,/reference,/dcsrsoft,/work,/scratch"

singularity run /dcsrsoft/singularity/containers/tetools-2.0.7.sif RepeatModeler -database TytoDB -threads 30 -srand 1 > TytoDB.out

echo ""
echo "End of the job [$(date)]"
echo ""
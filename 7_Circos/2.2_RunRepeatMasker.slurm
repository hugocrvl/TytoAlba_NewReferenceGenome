#!/bin/bash

#SBATCH --job-name=2.2_RunRepeatMasker.slurm
#SBATCH --output=2.2_RunRepeatMasker_%A.out

#SBATCH --mail-type=ALL
#SBATCH --mail-user=hugo.corval@unil.ch

#SBATCH --cpus-per-task=30
#SBATCH --mem=15G
#SBATCH --time=02:00:00

echo ""
echo "Start of the job [$(date)]"
echo ""

if [ $# -ne 3 ]; then
  echo "Usage: $0 FASTA REPEATS_FILE OUTPUT_DIR" >&2
  exit 1
fi

REF=${1}
REPEATS_FILE=${2}
OUTPUT_DIR=${3}

module load singularityce

export SINGULARITY_BINDPATH="/users,/reference,/dcsrsoft,/work,/scratch"

# singularity run /dcsrsoft/singularity/containers/tetools-2.0.7.sif RepeatMasker -pa 30 -small -xsmall -dir ${OUTPUT_DIR} -lib ${REPEATS_FILE} ${REF}
singularity run /dcsrsoft/singularity/containers/tetools-2.0.7.sif RepeatMasker -pa 30 -dir ${OUTPUT_DIR} -lib ${REPEATS_FILE} ${REF}

echo ""
echo "End of the job [$(date)]"
echo ""
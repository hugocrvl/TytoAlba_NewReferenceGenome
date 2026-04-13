#!/bin/bash

#SBATCH --job-name=7_MappingCDS2UNZIP.slurm
#SBATCH --output=7.1_MappingCDS2UNZIP_%A.out

#SBATCH --mail-type=ALL
#SBATCH --mail-user=hugo.corval@unil.ch

#SBATCH --cpus-per-task=5
#SBATCH --mem=20G
#SBATCH --time=00:30:00

if [ $# -ne 3 ]; then
  echo "Usage: $0 ASSEMBLY.fasta CDS.fasta OUTPUT.sam" >&2
  exit 1
fi

echo ""
echo "Start of the job [$(date)]"
echo ""

ASSEMBLY=${1}
CDS=${2}
OUTPUT=${3}

module load minimap2

echo "Mapping of CDS started [...]"

minimap2 -t 5 -x splice --secondary=no --cs=long ${ASSEMBLY} ${CDS} > ${OUTPUT}

echo "Done mapping"

echo ""
echo "End of the job: [$(date)]"
echo ""
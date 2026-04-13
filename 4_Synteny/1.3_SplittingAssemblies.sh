#!/bin/bash

# Creating the output directory
mkdir -p ${1}

# Setting variables
assembly=${2}
listLGs=${3}

# Splitting assembly per LG
module load samtools
for LG in $(cat ${listLGs}); do
    echo "${LG}"
    samtools faidx ${assembly} ${LG} > ${1}/${LG}.fasta
done
#!/usr/bin/env bash
# bash RunSnpable.sh <Genome.fa> <kmerSize> <OutputPrefix>

set -euo pipefail

date


module load bwa

#### Inputs ####
GENOME=$1
K=$2
PREFIX=$3

#### Paths to tools ####
SEQBILITY_DIR="/users/tcumer/WORK/TOOLS/seqbility-20091110"
SPLITFA="${SEQBILITY_DIR}/splitfa"
GEN_RAW_MASK="${SEQBILITY_DIR}/gen_raw_mask.pl"
GEN_MASK="${SEQBILITY_DIR}/gen_mask"

THREADS=4
WORKDIR="${PREFIX}_tmp"
mkdir -p "${WORKDIR}"

echo "Starting extraction of overlapping ${K}-mer subsequences from ${GENOME}"

# Stream all k-mers into a single FASTA file (no chunking needed unless you want parallel alignment)
"${SPLITFA}" "${GENOME}" "${K}" > "${PREFIX}_split.${K}.fa"

echo "Checking BWA index for ${GENOME}"
if [ -f "${GENOME}.bwt" ]; then
    echo "${GENOME} already indexed"
else
    echo "Indexing ${GENOME}"
    bwa index "${GENOME}"
fi

echo "Aligning ${K}-mer reads to the genome with BWA (aln/samse) and generating SAM"
bwa aln -t "${THREADS}" -l 32 -k 2 -O 3 -E 3 "${GENOME}" "${PREFIX}_split.${K}.fa" > "${PREFIX}_split.${K}.sai"
bwa samse "${GENOME}" "${PREFIX}_split.${K}.sai" "${PREFIX}_split.${K}.fa" > "${PREFIX}_split.${K}.sam"

echo "Reads aligned, generating rawMask"
"${GEN_RAW_MASK}" "${PREFIX}_split.${K}.sam" > "${PREFIX}_rawMask.${K}.fa"

echo "Raw mask created as ${PREFIX}_rawMask.${K}.fa, now generating final masks with different stringency"
"${GEN_MASK}" -l "${K}" -r 0.50 "${PREFIX}_rawMask.${K}.fa" > "${PREFIX}_mask.${K}.50.fa"
"${GEN_MASK}" -l "${K}" -r 0.90 "${PREFIX}_rawMask.${K}.fa" > "${PREFIX}_mask.${K}.90.fa"
"${GEN_MASK}" -l "${K}" -r 0.95 "${PREFIX}_rawMask.${K}.fa" > "${PREFIX}_mask.${K}.95.fa"
"${GEN_MASK}" -l "${K}" -r 0.99 "${PREFIX}_rawMask.${K}.fa" > "${PREFIX}_mask.${K}.99.fa"

echo "Cleaning up temporary files"
rm -f "${PREFIX}_split.${K}.sai" || true
# If you used chunking:
# rm -rf "${WORKDIR}"

echo "All done!"
date


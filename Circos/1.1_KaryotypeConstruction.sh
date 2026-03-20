#!/bin/bash
# Create a Circos karyotype file from a FASTA
# Usage:
#   ./fasta_to_karyotype.sh input.fasta COLOR > karyotype.txt
# Examples:
#   ./fasta_to_karyotype.sh genome.fa blue > karyotype.Father.txt
#   ./fasta_to_karyotype.sh contigs.fa vdgrey > karyotype.Mother.txt
#
# Notes:
# - Handles multi-line/wrapped FASTA sequences.
# - Uses the header (after '>') up to the first whitespace as the sequence name.
# - START is always 0, END is total sequence length in bases.
# - If input is .gz, the script will auto-decompress via zcat.
# - Sequences are expected to be named chr1, chr2, ..., chrZ.
set -euo pipefail
if [ $# -ne 2 ]; then
  echo "Usage: $0 <FASTA> <COLOR>" >&2
  exit 1
fi
FASTA="$1"
COLOR="$2"
if [ ! -f "$FASTA" ]; then
  echo "Error: FASTA file '$FASTA' not found." >&2
  exit 1
fi
READER="cat"
case "$FASTA" in
  *.gz) READER="zcat" ;;
esac
$READER "$FASTA" | awk -v color="$COLOR" '
  /^>/ {
    if (seqname != "" && seqlen > 0) {
      label = seqname
      sub(/^.*chr/, "", label)
      printf("chr - %s %s 0 %d %s\n", seqname, label, seqlen, color)
    }
    header = substr($0, 2)
    split(header, a, /[ \t]/)
    seqname = a[1]
    seqlen = 0
    next
  }
  {
    gsub(/[ \t\r]/, "", $0)
    if (length($0) > 0) {
      seqlen += length($0)
    }
  }
  END {
    if (seqname != "" && seqlen > 0) {
      label = seqname
      sub(/^.*chr/, "", label)
      printf("chr - %s %s 0 %d %s\n", seqname, label, seqlen, color)
    }
  }
'
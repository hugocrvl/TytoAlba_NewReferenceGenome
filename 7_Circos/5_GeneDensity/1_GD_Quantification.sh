#!/bin/bash
# Usage:
#   bash 5_GeneDensity/1_GD_Quantification.sh <genes.gff_or_txt> [WIN]
# Example:
#   bash 5_GeneDensity/1_GD_Quantification.sh genes.txt > GD_output.txt

set -euo pipefail

GFF="${1:-}"
WIN="${2:-100000}"

if [[ -z "$GFF" || ! -f "$GFF" ]]; then
  echo "ERROR: Input file not found: $GFF" >&2
  exit 1
fi

# Normalize CRLF to LF if dos2unix is installed
if command -v dos2unix >/dev/null 2>&1; then
  tmpfile=$(mktemp)
  cp "$GFF" "$tmpfile"
  dos2unix -q "$tmpfile" >/dev/null 2>&1 || true
  GFF="$tmpfile"
fi

awk -v W="$WIN" '
BEGIN {
  FS = "[ \t]+"; OFS = "\t"
}

# Skip header or commented lines
NR==1 && ($1=="seqid" || $1=="#seqid") { next }
$0 ~ /^#/ { next }

# Process only lines where type == "gene"
$3 != "gene" { next }

{
  chrom = $1
  start = $5 + 0
  end   = $6 + 0

  # Remove trailing CR (safe cleanup)
  sub(/\r$/, "", chrom)

  # Validate coordinates
  if (start <= 0 || end <= 0) next
  if (end < start) { tmp=start; start=end; end=tmp }

  # Window indexing
  bin_start = int((start - 1) / W)
  bin_end   = int((end   - 1) / W)

  for (b = bin_start; b <= bin_end; b++) {
    key = chrom ":" b
    counts[key]++
    bins[key] = chrom OFS (b*W + 1) OFS ((b+1)*W)
  }
}

END {
  for (key in bins) {
    print bins[key], counts[key]
  }
}
' "$GFF" | sort -k1,1 -k2,2n
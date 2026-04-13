
#!/usr/bin/env bash
# Compute mean of the 4th column in a per-base bedGraph over 100,000 bp windows.
# Output: tab-delimited text file with columns: chrom, win_start, win_end, mean
#
# Usage:
#   bedgraph_mean_100kbp.sh input.bedgraph > output.tsv

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 input.bedgraph" >&2
  exit 1
fi

IN="$1"
WIN=100000

awk -v W="$WIN" 'BEGIN{
  OFS="\t"
}
NF < 4 || $0 ~ /^#/ { next }

{
  chr = $1
  s   = $2 + 0
  e   = $3 + 0
  v   = $4 + 0

  if (e <= s) next

  win = int(s / W)
  key = chr SUBSEP win

  sum[key] += v
  count[key] += (e - s)  # per-base rows, so this is 1
}

END{
  PROCINFO["sorted_in"] = "@ind_str_asc"

  for (k in sum) {
    split(k, arr, SUBSEP)
    chr = arr[1]
    win = arr[2] + 0
    win_start = win * W
    win_end   = win_start + W

    mean_cov = sum[k] / count[k]  # mean over covered bases only
    print chr, win_start, win_end, mean_cov
  }
}
' "$IN"

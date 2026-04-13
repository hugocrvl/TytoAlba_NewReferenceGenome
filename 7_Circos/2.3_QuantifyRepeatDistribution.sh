
#!/usr/bin/env bash
set -euo pipefail

# rm_out_to_windows_100kb.sh
# From RepeatMasker .out, compute fraction of masked bases per 100 kb window.
#
# Usage:
#   ./rm_out_to_windows_100kb.sh RepeatMasker.out [genome.fai] > out.tsv
#
# Output (TSV, no header):
#   chrom  win_start  win_end  fraction_masked  (fraction in 0..1, e.g., 0.010 = 1%)
#
# Window size:
WIN=100000

if [ $# -lt 1 ] || [ $# -gt 2 ]; then
  echo "Usage: $0 RepeatMasker.out [genome.fai]" >&2
  exit 1
fi

OUT="$1"
FAI="${2:-}"

# Pipeline:
# 1) Parse .out and split intervals across 100 kb windows
# 2) Sort
# 3) Merge overlaps and sum covered bp per (chrom, window_start)
# 4) Compute window_end and fraction; load .fai if provided via getline (no brace-pipe)

awk -v WIN="$WIN" '
  BEGIN { OFS = "\t" }
  # Only process lines where first field is numeric (skip header lines)
  $1 ~ /^[0-9]+$/ {
    chrom = $5
    s = $6 + 0
    e = $7 + 0
    if (s > e) { t = s; s = e; e = t }
    if (s < 1) s = 1

    # Emit segments split by 100 kb windows
    while (s <= e) {
      wstart = int((s - 1) / WIN) * WIN + 1
      wend   = wstart + WIN - 1
      segend = (e < wend ? e : wend)
      print chrom, wstart, s, segend
      s = segend + 1
    }
  }
' "$OUT" \
| sort -k1,1 -k2,2n -k3,3n \
| awk '
  BEGIN { OFS = "\t" }
  {
    key = $1 FS $2
    if (key != prev_key) {
      if (prev_key != "") {
        covered += cur_end - cur_start + 1
        print prev_chr, prev_ws, covered
      }
      prev_key = key
      prev_chr = $1
      prev_ws  = $2 + 0
      covered  = 0
      cur_start = $3 + 0
      cur_end   = $4 + 0
    } else {
      if ($3 <= cur_end + 1) {
        if ($4 > cur_end) cur_end = $4 + 0
      } else {
        covered += cur_end - cur_start + 1
        cur_start = $3 + 0
        cur_end   = $4 + 0
      }
    }
  }
  END {
    if (prev_key != "") {
      covered += cur_end - cur_start + 1
      print prev_chr, prev_ws, covered
    }
  }
' \
| awk -v WIN="$WIN" -v FAI="$FAI" '
  BEGIN {
    OFS = "\t"
    has_fai = (FAI != "")
    if (has_fai) {
      # Load contig lengths from FAI: <name> <length> ...
      while ((getline line < FAI) > 0) {
        split(line, a, /\t/)
        len[a[1]] = a[2] + 0
      }
      close(FAI)
    }
  }
  {
    chrom = $1
    wstart = $2 + 0
    masked = $3 + 0
    wend = wstart + WIN - 1
    if (has_fai && (chrom in len) && wend > len[chrom]) {
      wend = len[chrom]
    }
    wlen = (has_fai && (chrom in len) ? (wend - wstart + 1) : WIN)
    frac = (wlen > 0 ? masked / wlen : 0.0)
    # 3 decimals like "0.010"
    printf "%s\t%d\t%d\t%.3f\n", chrom, wstart, wend, frac
  }
'

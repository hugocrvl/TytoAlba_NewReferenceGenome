#!/usr/bin/env bash
set -euo pipefail

# red_rpt_to_windows_100kb.sh
# From RED .rpt, compute fraction of masked bases per 100 kb window.
#
# Format expected in .rpt:
#   >chrom:start-end
# Example:
#   >Fa_LG01:0-57838
#
# Usage:
#   ./red_rpt_to_windows_100kb.sh repeats.rpt [genome.fai] [--all] > out.tsv
#
# Output (TSV, no header):
#   chrom  win_start  win_end  fraction_masked
#     - fraction in 0..1 (e.g., 0.010 = 1%)
#
# Notes:
#   - Windows are 1-based, inclusive: 1-100000, 100001-200000, ...
#   - With [--all], emit all windows (including those with 0 masked bp); requires genome.fai.
#   - Without [--all], only windows having masked bases are printed.

WIN=100000

if [ $# -lt 1 ] || [ $# -gt 3 ]; then
  echo "Usage: $0 repeats.rpt [genome.fai] [--all]" >&2
  exit 1
fi

RPT="$1"
FAI=""
ALL="0"

if [ $# -ge 2 ]; then
  if [ "${2:-}" = "--all" ]; then
    ALL="1"
  else
    FAI="$2"
  fi
fi

if [ $# -ge 3 ]; then
  if [ "${3:-}" = "--all" ]; then
    ALL="1"
  else
    echo "Unknown third argument: $3 (did you mean --all?)" >&2
    exit 1
  fi
fi

if [ "$ALL" = "1" ] && [ -z "$FAI" ]; then
  echo "[--all] requires genome.fai to know contig lengths." >&2
  exit 1
fi

# Pipeline:
# 1) Parse .rpt and split intervals across 100 kb windows
# 2) Sort
# 3) Merge overlaps and sum covered bp per (chrom, window_start)
# 4) Compute window_end and fraction; load .fai if provided via getline
# 5) Optionally (with --all) enumerate all windows using .fai and fill zeros

# Steps 1–3: produce per-window covered-bp totals
COVERED_STREAM="$(
  awk -v WIN="$WIN" '
    BEGIN { OFS = "\t" }
    # Lines like: >chrom:start-end
    $0 ~ /^>/ {
      line = $0
      sub(/^>/, "", line)
      n = split(line, a, /[:\-]/)  # a[1]=chrom, a[2]=start, a[3]=end
      if (n < 3) next
      chrom = a[1]
      s = a[2] + 0
      e = a[3] + 0

      # Normalize
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
  ' "$RPT" \
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
        # Merge within-window overlaps
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
  '
)"

if [ "$ALL" = "0" ]; then
  # Step 4 (no fill): compute fraction for windows with any coverage
  awk -v WIN="$WIN" -v FAI="$FAI" '
    BEGIN {
      OFS = "\t"
      has_fai = (FAI != "")
      if (has_fai) {
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
      printf "%s\t%d\t%d\t%.3f\n", chrom, wstart, wend, frac
    }
  ' <<< "$COVERED_STREAM"
else
  # Step 5 (fill all windows): requires FAI
  awk -v WIN="$WIN" -v FAI="$FAI" '
    BEGIN {
      OFS = "\t"
      # Load contig lengths
      while ((getline line < FAI) > 0) {
        split(line, a, /\t/)
        len[a[1]] = a[2] + 0
        order[++nchr] = a[1]  # preserve FAI order
      }
      close(FAI)
    }
    # Bring in covered-bp per (chrom, wstart)
    {
      cov[$1, $2] = $3 + 0
    }
    END {
      for (i = 1; i <= nchr; i++) {
        chr = order[i]
        L = len[chr]
        if (L <= 0) continue
        for (wstart = 1; wstart <= L; wstart += WIN) {
          wend = wstart + WIN - 1
          if (wend > L) wend = L
          wlen = wend - wstart + 1
          masked = ((chr, wstart) in cov ? cov[chr, wstart] : 0)
          frac = (wlen > 0 ? masked / wlen : 0.0)
          printf "%s\t%d\t%d\t%.3f\n", chr, wstart, wend, frac
        }
      }
    }
  ' <<< "$COVERED_STREAM" | sort -k1,1 -k2,2n
fi
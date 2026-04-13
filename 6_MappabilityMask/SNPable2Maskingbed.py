#!/usr/bin/env python3
import sys, gzip, re

def lines(fp):
    """Yield decoded text lines from file-like object."""
    for b in fp:
        yield b.decode() if isinstance(b, bytes) else b

def open_any(path):
    """Open stdin, plain text, or gzipped file transparently."""
    if path in ("-", None):
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "r")

def main(path):
    chrom = None
    pos = 0              # 0-based position counter
    in_run = False       # currently inside a run of NOT-3 bases?
    start = -1

    with open_any(path) as fh:
        for line in lines(fh):
            if line.startswith(">"):
                # flush ongoing run before switching chromosome
                if in_run:
                    print(chrom, start, pos, sep="\t")
                    in_run = False

                chrom = line[1:].strip().split()[0]
                pos = 0
                continue

            # normalize the line: remove whitespace
            seq = re.sub(r"\s+", "", line)
            if not seq:
                continue

            for ch in seq:
                # any character except "3" counts as "non-SNPable"
                if ch != "3":
                    if not in_run:
                        start = pos
                        in_run = True
                else:
                    if in_run:
                        print(chrom, start, pos, sep="\t")
                        in_run = False
                pos += 1

        # flush final run at end of file
        if in_run:
            print(chrom, start, pos, sep="\t")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "-")

#!/usr/bin/env python3
import sys, gzip, re

def lines(fp):
    for b in fp:
        yield b.decode() if isinstance(b, bytes) else b

def open_any(path):
    if path in ('-', None):
        return sys.stdin
    if path.endswith('.gz'):
        return gzip.open(path, 'rb')
    return open(path, 'r')

def main(path):
    chrom = None
    pos = 0
    in_run = False
    start = -1

    with open_any(path) as fh:
        for line in lines(fh):
            if line.startswith('>'):
                if in_run:
                    print(chrom, start, pos, sep='\t')
                    in_run = False
                chrom = line[1:].strip().split()[0]
                pos = 0
                continue
            seq = re.sub(r'\s+', '', line)
            for ch in seq:
                if ch == '3':
                    if not in_run:
                        start = pos; in_run = True
                else:
                    if in_run:
                        print(chrom, start, pos, sep='\t')
                        in_run = False
                pos += 1
        if in_run:
            print(chrom, start, pos, sep='\t')

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '-')

#!/usr/bin/perl
#
# random list of selected reads (argument 1) will be chopped in kmers
# only kmers without N will be printed
# they are printed in lexically sorted orientation (e.g. revcomp if needed)
#
# (c) N.Guex 2024-2025
# Bioinformatics Competence Center (BICC)
# University of Lausanne
# Switzerland
# https://bix.unil.ch
#
# NG 20240213

use strict;
my $kmerLen = $ARGV[0];
my $fn = $ARGV[1];
my @desiredLine;
my $l = 0;

# ----- load lines to consider
open F, "$fn" or die "cannot open $fn\n";
while(<F>)
{
    chomp;
    $desiredLine[$l++] = $_;
}
close F;

# ----------- print kmers of desired lines
my $line = 1;
my $idx = 0;
my $validKmersCnt = 0;
while(<STDIN>)
{
    my $SEQ = <STDIN>;
    <STDIN>;
    <STDIN>; # quality
    if ($line == $desiredLine[$idx])
    {
        chomp($SEQ);
        for (my $i=0; $i<= (length($SEQ)-$kmerLen); $i++) {
            my $kmer = substr($SEQ,$i,$kmerLen);
            unless ($kmer =~ /N/)
      {
        my $REV = reverse($kmer);
        $REV =~ tr/ACGTacgt/TGCAtgca/;
        if ($kmer le $REV) { print ">$line\n$kmer\n"; } else { print ">$line\n$REV\n"; }
        $validKmersCnt++;
      }
        }
        $idx++;
        if ($idx == $l) { last; }
    }
    $line++;
}
#-----------------------
print STDERR "$fn:$validKmersCnt\n";


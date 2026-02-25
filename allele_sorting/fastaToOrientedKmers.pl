#!/usr/bin/perl
# This takes a fasta file (sequence must be on one line) 
# It generates all kmers (only kmers without N will be printed)
# they are printed in lexically sorted orientation (e.g. revcomp if needed)
#
# (c) N.Guex 2024-2025
# Bioinformatics Competence Center (BICC)
# University of Lausanne
# Switzerland
# https://bix.unil.ch
#
# NG 20240213

$DBLEN=$ARGV[0];
$STEP=$ARGV[1];
$HDR=<STDIN>;
$SEQ=<STDIN>;
$SEQ =~ s/\-//g;
$len = length($SEQ);
if ($STEP eq "") { die "usage DBLEN STEP\n"; }
for ($i=0; $i<($len-$DBLEN); $i+=$STEP)
{
	$PAT = substr($SEQ,$i,$DBLEN);
        $REV = reverse($PAT);
	$REV =~ tr/ACGTacgt/TGCAtgca/;
	if ($PAT le $REV)
	{
		print "$PAT\n";
	}
       	else 
	{ 
		print "$REV\n";
	}
}


#!/usr/bin/perl
#
# This loads kmer counts specific for one parent and records the number
# of occurrences, and their total count to allow subsequent vote.
# expects Father counts in first column and Mother counts in second column
# # NG 20240213
#

use strict;
my $Fscore = 0;
my $Mscore = 0;
my $Fcnt = 0;
my $Mcnt = 0;

while(<STDIN>)
{
    if ($_ eq "0\t0\n") { next; }
    chomp;
    my (@X) = split("\t",$_);
    if ($X[1] == 0)
    {
        $Fscore += $X[0];
        $Fcnt++;
    }
    else
    {
        $Mscore += $X[1];
        $Mcnt++;
    }
}
print "$Fcnt,$Mcnt,$Fscore,$Mscore\n";


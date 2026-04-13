# This scipt collects KMERS present in parental illumina libraries
# kmers that are fed to jellyfish are already filtered and lexically oriented.
# it is important to NOT use the option -C when building the database for consistency during subsequent merging and handling of databases.
#
# (c) N.Guex 2024-2025
# Bioinformatics Competence Center (BICC)
# University of Lausanne
# Switzerland
# https://bix.unil.ch
#
# NG 20240213

KMER=49

time  ./jellyfish count -o  mother.jf${KMER}.L1.Uinf.selection1.bin  -m $KMER -c 4 -t 24 -L1 -s 8589934592 <(gzip -dc /tmp/nguex/Goudet/M028977_R1_paired.fq.gz        | ./collectKmerWithoutN.pl $KMER mother_selection1.txt ; gzip -dc /tmp/nguex/Goudet/M028977_R2_paired.fq.gz        | ../collectKmerWithoutN.pl $KMER mother_selection1.txt)
time  ./jellyfish count -o  mother.jf${KMER}.L1.Uinf.selection2.bin  -m $KMER -c 4 -t 24 -L1 -s 8589934592 <(gzip -dc /tmp/nguex/Goudet/M028977_R1_paired.fq.gz        | ./collectKmerWithoutN.pl $KMER mother_selection2.txt ; gzip -dc /tmp/nguex/Goudet/M028977_R2_paired.fq.gz        | ../collectKmerWithoutN.pl $KMER mother_selection2.txt)
time  ./jellyfish count -o  father.jf${KMER}.L1.Uinf.selection1.bin  -m $KMER -c 4 -t 24 -L1 -s 8589934592 <(gzip -dc /tmp/nguex/Goudet/M032330H_40_L?_R1_paired.fq.gz | ./collectKmerWithoutN.pl $KMER father_selection1.txt ; gzip -dc /tmp/nguex/Goudet/M032330H_40_L?_R2_paired.fq.gz | ../collectKmerWithoutN.pl $KMER father_selection1.txt)
time  ./jellyfish count -o  father.jf${KMER}.L1.Uinf.selection2.bin  -m $KMER -c 4 -t 24 -L1 -s 8589934592 <(gzip -dc /tmp/nguex/Goudet/M032330H_40_L?_R1_paired.fq.gz | ./collectKmerWithoutN.pl $KMER father_selection2.txt ; gzip -dc /tmp/nguex/Goudet/M032330H_40_L?_R2_paired.fq.gz | ../collectKmerWithoutN.pl $KMER father_selection2.txt)
time  ./jellyfish count -o  father.jf${KMER}.L1.Uinf.selection3.bin  -m $KMER -c 4 -t 24 -L1 -s 8589934592 <(gzip -dc /tmp/nguex/Goudet/M032330H_40_L?_R1_paired.fq.gz | ./collectKmerWithoutN.pl $KMER father_selection3.txt ; gzip -dc /tmp/nguex/Goudet/M032330H_40_L?_R2_paired.fq.gz | ../collectKmerWithoutN.pl $KMER father_selection3.txt)

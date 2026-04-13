# This scipt prepare databases with kmers present in one parent but NOT in the other using the modified jellyfish version
# To balance the uneven sizes of databases, consider always only "two parent databases for the exclusion".
#
# (c) N.Guex 2024-2025
# Bioinformatics Competence Center (BICC)
# University of Lausanne
# Switzerland
# https://bix.unil.ch
#
# NG 20240213

KMER=49
LOWEST_IN=L1
LOWEST_OUT=L2

# make Mother only databases
time  ./bin/jellyfish merge -o only_mother.jf${KMER}.${LOWEST_OUT}.Uinf.selection1.NOTfather1and2.bin  -${LOWEST_OUT} -O  mother.jf${KMER}.${LOWEST_IN}.Uinf.selection1.bin          father.jf${KMER}.${LOWEST_IN}.Uinf.selection1.bin father.jf${KMER}.${LOWEST_IN}.Uinf.selection2.bin
time  ./bin/jellyfish merge -o only_mother.jf${KMER}.${LOWEST_OUT}.Uinf.selection2.NOTfather1and2.bin  -${LOWEST_OUT} -O  mother.jf${KMER}.${LOWEST_IN}.Uinf.selection2.bin          father.jf${KMER}.${LOWEST_IN}.Uinf.selection1.bin father.jf${KMER}.${LOWEST_IN}.Uinf.selection2.bin

time  ./bin/jellyfish merge -o only_mother.jf${KMER}.${LOWEST_OUT}.Uinf.selection1.NOTfather1and3.bin  -${LOWEST_OUT} -O  mother.jf${KMER}.${LOWEST_IN}.Uinf.selection1.bin          father.jf${KMER}.${LOWEST_IN}.Uinf.selection1.bin father.jf${KMER}.${LOWEST_IN}.Uinf.selection3.bin
time  ./bin/jellyfish merge -o only_mother.jf${KMER}.${LOWEST_OUT}.Uinf.selection2.NOTfather1and3.bin  -${LOWEST_OUT} -O  mother.jf${KMER}.${LOWEST_IN}.Uinf.selection2.bin          father.jf${KMER}.${LOWEST_IN}.Uinf.selection1.bin father.jf${KMER}.${LOWEST_IN}.Uinf.selection3.bin

time  ./bin/jellyfish merge -o only_mother.jf${KMER}.${LOWEST_OUT}.Uinf.selection1.NOTfather2and3.bin  -${LOWEST_OUT} -O  mother.jf${KMER}.${LOWEST_IN}.Uinf.selection1.bin          father.jf${KMER}.${LOWEST_IN}.Uinf.selection2.bin father.jf${KMER}.${LOWEST_IN}.Uinf.selection3.bin
time  ./bin/jellyfish merge -o only_mother.jf${KMER}.${LOWEST_OUT}.Uinf.selection2.NOTfather2and3.bin  -${LOWEST_OUT} -O  mother.jf${KMER}.${LOWEST_IN}.Uinf.selection2.bin          father.jf${KMER}.${LOWEST_IN}.Uinf.selection2.bin father.jf${KMER}.${LOWEST_IN}.Uinf.selection3.bin



# make Father only databases
time  ./bin/jellyfish merge -o only_father.jf${KMER}.${LOWEST_OUT}.Uinf.selection1.NOTmother1and2.bin  -${LOWEST_OUT} -O  father.jf${KMER}.${LOWEST_IN}.Uinf.selection1.bin          mother.jf${KMER}.${LOWEST_IN}.Uinf.selection1.bin mother.jf${KMER}.${LOWEST_IN}.Uinf.selection2.bin
time  ./bin/jellyfish merge -o only_father.jf${KMER}.${LOWEST_OUT}.Uinf.selection2.NOTmother1and2.bin  -${LOWEST_OUT} -O  father.jf${KMER}.${LOWEST_IN}.Uinf.selection2.bin          mother.jf${KMER}.${LOWEST_IN}.Uinf.selection1.bin mother.jf${KMER}.${LOWEST_IN}.Uinf.selection2.bin
time  ./bin/jellyfish merge -o only_father.jf${KMER}.${LOWEST_OUT}.Uinf.selection3.NOTmother1and2.bin  -${LOWEST_OUT} -O  father.jf${KMER}.${LOWEST_IN}.Uinf.selection3.bin          mother.jf${KMER}.${LOWEST_IN}.Uinf.selection1.bin mother.jf${KMER}.${LOWEST_IN}.Uinf.selection2.bin


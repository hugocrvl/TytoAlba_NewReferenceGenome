# This scipt will make all possible parental specific kmer database comparisons and gather the corresponding counts for each read
#
# (c) N.Guex 2024-2025
# Bioinformatics Competence Center (BICC)
# University of Lausanne
# Switzerland
# https://bix.unil.ch
#
# NG 20250707
#
#  cd /home/nguex/GOUDET/allele_sort/jf49o
#  cp -p *NOT*.bin /dev/shm/
#  split -l 23501 <(cut -f1 /data2/chris/tyto/ON_25k_table.txt)  # cut in 48 jobs
#  ln -s predict_ON_ReadParent25k.sh shmQueryL1.sh
#  for i in `ls -1 x??.sh`; do ./$i & done
#

JFBIN=./bin/jellyfish
DBDIR=/dev/shm
KMER=49
LOWEST_OUT=L2

DB=`echo $1 | cut -c1-3`
samtools faidx -i /data2/chris/tyto/${DB}_25k.fq.gz $1 |  ./fastaToOrientedKmers.pl $KMER 7      > /dev/shm/jf7.$1.00
cat /dev/shm/jf7.$1.00 | $JFBIN query -i $DBDIR/only_father.jf${KMER}.${LOWEST_OUT}.Uinf.selection1.NOTmother1and2.bin > /dev/shm/jf7.$1.F1notM12
cat /dev/shm/jf7.$1.00 | $JFBIN query -i $DBDIR/only_father.jf${KMER}.${LOWEST_OUT}.Uinf.selection2.NOTmother1and2.bin > /dev/shm/jf7.$1.F2notM12
cat /dev/shm/jf7.$1.00 | $JFBIN query -i $DBDIR/only_father.jf${KMER}.${LOWEST_OUT}.Uinf.selection3.NOTmother1and2.bin > /dev/shm/jf7.$1.F3notM12


cat /dev/shm/jf7.$1.00 | $JFBIN query -i $DBDIR/only_mother.jf${KMER}.${LOWEST_OUT}.Uinf.selection1.NOTfather1and2.bin > /dev/shm/jf7.$1.M1notF12
cat /dev/shm/jf7.$1.00 | $JFBIN query -i $DBDIR/only_mother.jf${KMER}.${LOWEST_OUT}.Uinf.selection2.NOTfather1and2.bin > /dev/shm/jf7.$1.M2notF12

cat /dev/shm/jf7.$1.00 | $JFBIN query -i $DBDIR/only_mother.jf${KMER}.${LOWEST_OUT}.Uinf.selection1.NOTfather1and3.bin > /dev/shm/jf7.$1.M1notF13
cat /dev/shm/jf7.$1.00 | $JFBIN query -i $DBDIR/only_mother.jf${KMER}.${LOWEST_OUT}.Uinf.selection2.NOTfather1and3.bin > /dev/shm/jf7.$1.M2notF13

cat /dev/shm/jf7.$1.00 | $JFBIN query -i $DBDIR/only_mother.jf${KMER}.${LOWEST_OUT}.Uinf.selection1.NOTfather2and3.bin > /dev/shm/jf7.$1.M1notF23
cat /dev/shm/jf7.$1.00 | $JFBIN query -i $DBDIR/only_mother.jf${KMER}.${LOWEST_OUT}.Uinf.selection2.NOTfather2and3.bin > /dev/shm/jf7.$1.M2notF23


# ----- GATHER COUNTS

# father1,2 or 3   compared to mother 1
compF01=`paste /dev/shm/jf7.$1.F1notM12 /dev/shm/jf7.$1.M1notF12 | ./chooseParent_F_M.pl`
compF02=`paste /dev/shm/jf7.$1.F2notM12 /dev/shm/jf7.$1.M1notF12 | ./chooseParent_F_M.pl`
compF03=`paste /dev/shm/jf7.$1.F3notM12 /dev/shm/jf7.$1.M1notF12 | ./chooseParent_F_M.pl`

# father1,2 or 3   compared to mother 2
compF04=`paste /dev/shm/jf7.$1.F1notM12 /dev/shm/jf7.$1.M2notF12 | ./chooseParent_F_M.pl`
compF05=`paste /dev/shm/jf7.$1.F2notM12 /dev/shm/jf7.$1.M2notF12 | ./chooseParent_F_M.pl`
compF06=`paste /dev/shm/jf7.$1.F3notM12 /dev/shm/jf7.$1.M2notF12 | ./chooseParent_F_M.pl`


# father1,2 or 3   compared to mother 1
compF07=`paste /dev/shm/jf7.$1.F1notM12 /dev/shm/jf7.$1.M1notF13 | ./chooseParent_F_M.pl`
compF08=`paste /dev/shm/jf7.$1.F2notM12 /dev/shm/jf7.$1.M1notF13 | ./chooseParent_F_M.pl`
compF09=`paste /dev/shm/jf7.$1.F3notM12 /dev/shm/jf7.$1.M1notF13 | ./chooseParent_F_M.pl`

# father1,2 or 3   compared to mother 2
compF10=`paste /dev/shm/jf7.$1.F1notM12 /dev/shm/jf7.$1.M2notF13 | ./chooseParent_F_M.pl`
compF11=`paste /dev/shm/jf7.$1.F2notM12 /dev/shm/jf7.$1.M2notF13 | ./chooseParent_F_M.pl`
compF12=`paste /dev/shm/jf7.$1.F3notM12 /dev/shm/jf7.$1.M2notF13 | ./chooseParent_F_M.pl`


# father1,2 or 3   compared to mother 1
compF13=`paste /dev/shm/jf7.$1.F1notM12 /dev/shm/jf7.$1.M1notF23 | ./chooseParent_F_M.pl`
compF14=`paste /dev/shm/jf7.$1.F2notM12 /dev/shm/jf7.$1.M1notF23 | ./chooseParent_F_M.pl`
compF15=`paste /dev/shm/jf7.$1.F3notM12 /dev/shm/jf7.$1.M1notF23 | ./chooseParent_F_M.pl`

# father1,2 or 3   compared to mother 2
compF16=`paste /dev/shm/jf7.$1.F1notM12 /dev/shm/jf7.$1.M2notF23 | ./chooseParent_F_M.pl`
compF17=`paste /dev/shm/jf7.$1.F2notM12 /dev/shm/jf7.$1.M2notF23 | ./chooseParent_F_M.pl`
compF18=`paste /dev/shm/jf7.$1.F3notM12 /dev/shm/jf7.$1.M2notF23 | ./chooseParent_F_M.pl`


#--------- RESULT

echo "#$1, $compF01, $compF02, $compF03, $compF04, $compF05, $compF06, $compF07, $compF08, $compF09, $compF10, $compF11, $compF12, $compF13, $compF14, $compF15, $compF16, $compF17, $compF18"

rm /dev/shm/jf7.$1.*


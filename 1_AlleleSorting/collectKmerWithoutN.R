# This scipt generate indices to randomly draw about the same number of reads from the illumina mother and father libraries
# It will generate 2 sets for mother and 3 sets for father
# These lists are used by the script collectKmerWithoutN.pl (called from collect.sh)
#
# (c) N.Guex 2024-2025
# Bioinformatics Competence Center (BICC)
# University of Lausanne
# Switzerland
# https://bix.unil.ch
#
# NG 20240213

set.seed(1)
s <- sample.int(93360939, size=46680469, replace = FALSE)
mother_selection1 <- sort(s)
mother_selection2 <- which(!(1:93360939 %in% mother_selection1))
write.table(file='mother_selection1.txt',mother_selection1,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(file='mother_selection2.txt',mother_selection2,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

s <- sample.int(172731951, size=2*46680469, replace = FALSE)
father_selection1 <- sort(s[1:46680469])
father_selection2 <- sort(s[-(1:46680469)])
father_selection3 <- which(!(1:172731951 %in% c(father_selection1,father_selection2)))[1:46680469]
write.table(file='father_selection1.txt',father_selection1,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(file='father_selection2.txt',father_selection2,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(file='father_selection3.txt',father_selection3,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)


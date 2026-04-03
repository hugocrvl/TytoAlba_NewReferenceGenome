# the physical_orders_full.txt file used here comes from topaloudis et al. 2024
# https://academic.oup.com/genetics/article/229/1/iyae190/7900905


## extract bed coordinates for each SNPs :
awk '!/#/ && NR>1 {print $2"\t"$3"\t"$3+200"\t"$2"_"$3"_"$1"_"$4}' 0_DATA/physical_orders_full.txt > 0_DATA/WindowSNPs.bed

## Extract the pseudo read for each SNPs in the previous assembly : 
ref='GCF_018691265.1_T.alba_DEE_v4.0_genomic.fna' # available at https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/691/265/GCF_018691265.1_T.alba_DEE_v4.0/

bedtools getfasta -fi ${ref} -bed 0_DATA/WindowSNPs.bed -fo 0_DATA/WindowSNPs.fasta -nameOnly

## Align it on the new refs : 
# Final_FatherAssembly.fa
# Final_MotherAssembly.fa

mkdir -p 1_AlignLG/


module load bwa samtools

NewRef_Fa="Final_FatherAssembly.fa"
bwa mem -M -t 4 $NewRef_Fa 0_DATA/WindowSNPs.fasta | samtools sort -@ 4 -O bam -o 1_AlignLG/LGSNPs_Fa.bam -
samtools view 1_AlignLG/LGSNPs_Fa.bam | awk '$5==60 {print $1,$3, $4}' > 1_AlignLG/LGSNPs_Fa.txt

# postprocess the file to generate the new linage maps
echo "#chr pos LG m" > Final_FatherAssembly_recombinationMap.txt
awk '{print $2" "$3}' 1_AlignLG/LGSNPs_Fa.txt > $$.1
awk -F "_" '{print $4" "$5}' 1_AlignLG/LGSNPs_Fa.txt > $$.2
paste $$.1 $$.2 | awk '{print $1, $2, $3, $4}' >> Final_FatherAssembly_recombinationMap.txt
rm $$.1 $$.2



NewRef_Mo="Final_MotherAssembly.fa"
bwa mem -M -t 4 $NewRef_Mo 0_DATA/WindowSNPs.fasta | samtools sort -@ 4 -O bam -o 1_AlignLG/LGSNPs_Mo.bam -
samtools view 1_AlignLG/LGSNPs_Mo.bam | awk '$5==60 {print $1,$3, $4}' > 1_AlignLG/LGSNPs_Mo.txt


# postprocess the file to generate the new linage maps
echo "#chr pos LG m" > Final_MotherAssembly_recombinationMap.txt
awk '{print $2" "$3}' 1_AlignLG/LGSNPs_Mo.txt > $$.1
awk -F "_" '{print $4" "$5}' 1_AlignLG/LGSNPs_Mo.txt > $$.2
paste $$.1 $$.2 | awk '{print $1, $2, $3, $4}' >> Final_MotherAssembly_recombinationMap.txt
rm $$.1 $$.2

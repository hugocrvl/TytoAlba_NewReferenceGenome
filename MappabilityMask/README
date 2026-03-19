# Prepare the reference files 

bwa index 0_Assemblies/Final_FatherAssembly.fa
samtools faidx 0_Assemblies/Final_FatherAssembly.fa

bwa index 0_Assemblies/Final_MotherAssembly.fa
samtools faidx 0_Assemblies/Final_MotherAssembly.fa

# Run SNPable pipeline

mkdir -p SNPable_OUT
bash RunSnpable.sh 0_Assemblies/Final_FatherAssembly.fa 150 SNPable_OUT/Snpable_Fa
bash RunSnpable.sh 0_Assemblies/Final_MotherAssembly.fa 150 SNPable_OUT/Snpable_Mo


# Convert the output into bed files

python SNPable2bed.py SNPable_OUT/Snpable_Fa_mask.150.50.fa > SNPable_OUT/Final_FatherAssembly.150.50.SNPableRegions.bed
python SNPable2Maskingbed.py SNPable_OUT/Snpable_Fa_mask.150.50.fa > SNPable_OUT/Final_FatherAssembly.150.50.MaskedRegions.bed

python SNPable2bed.py SNPable_OUT/Snpable_Fa_mask.150.90.fa > SNPable_OUT/Final_FatherAssembly.150.90.SNPableRegions.bed
python SNPable2Maskingbed.py SNPable_OUT/Snpable_Fa_mask.150.90.fa > SNPable_OUT/Final_FatherAssembly.150.90.mask.bed

python SNPable2bed.py SNPable_OUT/Snpable_Fa_mask.150.95.fa > SNPable_OUT/Final_FatherAssembly.150.95.SNPableRegions.bed
python SNPable2Maskingbed.py SNPable_OUT/Snpable_Fa_mask.150.95.fa > SNPable_OUT/Final_FatherAssembly.150.95.mask.bed

python SNPable2bed.py SNPable_OUT/Snpable_Fa_mask.150.99.fa > SNPable_OUT/Final_FatherAssembly.150.99.SNPableRegions.bed
python SNPable2Maskingbed.py SNPable_OUT/Snpable_Fa_mask.150.99.fa > SNPable_OUT/Final_FatherAssembly.150.99.mask.bed


python SNPable2bed.py SNPable_OUT/Snpable_Mo_mask.150.50.fa > SNPable_OUT/Final_MotherAssembly.150.50.SNPableRegions.bed
python SNPable2Maskingbed.py SNPable_OUT/Snpable_Mo_mask.150.50.fa > SNPable_OUT/Final_MotherAssembly.150.50.mask.bed

python SNPable2bed.py SNPable_OUT/Snpable_Mo_mask.150.90.fa > SNPable_OUT/Final_MotherAssembly.150.90.SNPableRegions.bed
python SNPable2Maskingbed.py SNPable_OUT/Snpable_Mo_mask.150.90.fa > SNPable_OUT/Final_MotherAssembly.150.90.mask.bed

python SNPable2bed.py SNPable_OUT/Snpable_Mo_mask.150.95.fa > SNPable_OUT/Final_MotherAssembly.150.95.SNPableRegions.bed
python SNPable2Maskingbed.py SNPable_OUT/Snpable_Mo_mask.150.95.fa > SNPable_OUT/Final_MotherAssembly.150.95.mask.bed

python SNPable2bed.py SNPable_OUT/Snpable_Mo_mask.150.99.fa > SNPable_OUT/Final_MotherAssembly.150.99.SNPableRegions.bed
python SNPable2Maskingbed.py SNPable_OUT/Snpable_Mo_mask.150.99.fa > SNPable_OUT/Final_MotherAssembly.150.99.mask.bed



# Categorisation of the chromosomes 

## Description

For the manuscript, we need to categorise chromosomes into three categories: macro-, micro- and dot-. Instead of only using the size as a criteria, we want to apply a PCA on several characteristics: GC, repeat, gene, size.

## Data gathering

GC- and repeats proportion and exonic length data have been computed in the `../Figures/Figure_4/DotChromosomes.R` script.

``` R
# I previously computed all the *.exons, *.repeats and GC.* data.frame using the dedicated script
index.father = read.table("../../0_FinalFiles/Final_FatherAssembly.fa.fai", header = F)
chromosomes.Father = data.frame(chr = index.father$V1,
                                size = Father.exons$chromosome.size, 
                                exons = Father.exons$relative.exonic.length,
                                repeats = Father.repeats$relative.repeats.length,
                                GC = GC.Father$GC)
write.table(chromosomes.Father, "../../15_ChromosomeCategories/ProportionsFather.txt", quote = F, row.names = F, col.names = T)

index.mother = read.table("../../0_FinalFiles/Final_MotherAssembly.fa.fai", header = F)
chromosomes.Mother = data.frame(chr = index.mother$V1,
                                size = Mother.exons$chromosome.size, 
                                exons = Mother.exons$relative.exonic.length,
                                repeats = Mother.repeats$relative.repeats.length,
                                GC = GC.Mother$GC)
write.table(chromosomes.Mother, "../../15_ChromosomeCategories/ProportionsMother.txt", quote = F, row.names = F, col.names = T)
```

## Genes features and repeat types per chromosome

Tristan did it, the results can be found here: `/work/FAC/FBM/DEE/jgoudet/barn_owl/tcumer/Analyses_Genome/D_NewRefGenome_22102025/8_ChromosomeFeatures/`.

``` bash
# Genes features
cp /work/FAC/FBM/DEE/jgoudet/barn_owl/tcumer/Analyses_Genome/D_NewRefGenome_22102025/8_ChromosomeFeatures/Father.complete.genomic.gff.stats ./GenesFeaturesFather.txt
cp /work/FAC/FBM/DEE/jgoudet/barn_owl/tcumer/Analyses_Genome/D_NewRefGenome_22102025/8_ChromosomeFeatures/Mother.complete.genomic.gff.stats ./GenesFeaturesMother.txt
# Repeats types
cp /work/FAC/FBM/DEE/jgoudet/barn_owl/tcumer/Analyses_Genome/D_NewRefGenome_22102025/8_ChromosomeFeatures/Final_FatherAssembly.fa.out.perChr ./RepeatsTypesFather.txt
cp /work/FAC/FBM/DEE/jgoudet/barn_owl/tcumer/Analyses_Genome/D_NewRefGenome_22102025/8_ChromosomeFeatures/Final_MotherAssembly.fa.out.perChr ./RepeatsTypesMother.txt
```

## Merging and creating the final data.frame

I used the `./2_CategorisationPCA.R` script.

## PCA

To perform the rest of the analysis, I wrote an R script, `./2_CategorisationPCA.R`.

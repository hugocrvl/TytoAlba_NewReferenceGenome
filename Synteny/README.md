
# Evaluation of the synteny between the paternal and maternal alleles

## Table of contents

1. [Description](#description)
2. [Input creatin](#creation-of-the-inputs)
3. [Alignment whole genome](#alignment)

## Description

I want to evaluate the synteny between parental haplotypes. I need to align each Mother LG to Father LG separately and then plot it using `SVbyEye`.

## Creation of the inputs

First, I need to split the assemblies per LG using samtools.

``` bash
# For the Father assembly
grep "^>" ../0_FinalFiles/Final_FatherAssembly.fa | sed 's/>//' > 1.1_Chromosomes_Father.txt
bash 1.3_SplittingAssemblies.sh 1.4_SplitAssembliesPerLGs_Father ../0_FinalFiles/Final_FatherAssembly.fa 1.1_Chromosomes_Father.txt

# For the Mother assembly
grep "^>" ../0_FinalFiles/Final_MotherAssembly.fa | sed 's/>//' > 1.2_Chromosomes_Mother.txt
bash 1.3_SplittingAssemblies.sh 1.5_SplitAssembliesPerLGs_Mother ../0_FinalFiles/Final_MotherAssembly.fa 1.2_Chromosomes_Mother.txt
```

## Alignment

I use `minimap2` to align every pair of parental LGs.

``` bash
sbatch 2_Alignment.slurm
```

I also aligned the two haplotypes entirely, without splitting per LG

``` bash
sbatch 3_AlignmentWholeGenome.slurm
```

Figure were created on my local computer using `R`.

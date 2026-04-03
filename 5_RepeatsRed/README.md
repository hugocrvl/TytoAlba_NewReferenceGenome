
# Detection of repeat regions across the new assembly

## Table of contents

1. [Description](#description)
2. [Red tool](#detection-with-red)
3. [Summarising information](#summarise-information-per-windows)

## Description

In this directory, I am going to investigate any toggling within the LG11 inversion. First, I quantify the amount of repeats across the assembly haplotypes with red.

## Detection with red

Tristan previously installed this tool on his directory.

``` bash
mkdir -p ./0.1_FaReference/
cp ../0_FinalFiles/Final_FatherAssembly.fa ./0.1_FaReference/
mkdir -p ./0.2_MoReference/
cp ../0_FinalFiles/Final_MotherAssembly.fa ./0.2_MoReference/

sbatch 1_Red.slurm
```

## Summarise information per windows

I want to quantify the percentage of repeats per windows of 100kbp. To do so, I have the `2_ProcessRedOutput.sh` script.

``` bash
bash 2_ProcessRedOutput.sh 1.3.1_Red_Repeats_Output_Fa/Final_FatherAssembly.rpt > 2.1_RepeatsFa_W100kbp.tsv
bash 2_ProcessRedOutput.sh 1.3.2_Red_Repeats_Output_Mo/Final_MotherAssembly.rpt > 2.2_RepeatsMo_W100kbp.tsv
```

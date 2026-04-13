
# Quality assessment of the reference genome

## Table of contents

- [Description](#description)
- [BlobTools2 installation](#blobtools2-installation)
- [FASTA files](#fasta-files)
- [NCBI input files](#ncbi-input-files)
- [BUSCO](#busco-installation-and-busco-score)
- [BlobTools2 pipeline](#creation-of-a-blobtools2-directory)
- [Snail plot](#snail-plot)
- [UNZIP](#unzip)
- [Quast](#quast)

## Description

In this directory, I will proceed to the quality assessment of our genome assembly. The goal would be to have stand-alone working script, with several input files that the user could change, for instance: `.fasta` file, `busco.tsv`, etc.

## BlobTools2 installation

Use of conda to install the package.

``` bash
module load miniforge3
conda_init
conda create -y -n btk -c conda-forge python=3.9
conda activate btk
pip install "blobtoolkit[full]"
blobtools -h
conda deactivate
```

## FASTA files

WARNING: I'm not using these FASTA files anymore, but the ones located in `/work/FAC/FBM/DEE/jgoudet/barn_owl/hcorval/Chapter_NewReferenceGenome/0_FinalFiles`.

Creation of a directory that contains the last version of the assembly we have.

``` bash

mkdir -p 1_FASTA
# Copying the father haplotype
cp /nas/FAC/FBM/DEE/jgoudet/barn_owl/D2c/bix/202512_final/all_Fa.fa ./1_FASTA
# Copying the mother haplotype
cp /nas/FAC/FBM/DEE/jgoudet/barn_owl/D2c/bix/202512_final/all_Mo.fa ./1_FASTA

# Indexing both fasta files
module load samtools
samtools faidx ./1_FASTA/all_Fa.fa
samtools faidx ./1_FASTA/all_Mo.fa
```

## NCBI input files

I download the NCBI dump containing the entire NCBI taxonomy database. At this point, I'm unsure whether it is useful or not.

``` bash
mkdir -p 2_NCBI
cd 2_NCBI
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
tar -xvzf new_taxdump.tar.gz
rm new_taxdump.tar.gz
```

For _Tyto alba_, I should use the `--taxid 56313` parameter.

## BUSCO installation and BUSCO score

I need to install `BUSCO` to have a BUSCO score for each assembly.

``` bash
module load miniforge3
conda_init
conda create -n busco
conda activate busco
conda install bioconda::busco
busco --help
conda deactivate
```

To run BUSCO on the reference genome:

``` bash
mkdir -p 3_BUSCO/
cd 3_BUSCO/
screen -S busco
# Father
bash 1_Run_BUSCO.sh ../../0_FinalFiles/Final_FatherAssembly.fa
# Mother
bash 1_Run_BUSCO.sh ../../0_FinalFiles/Final_MotherAssembly.fa
```

## Creation of a BlobTools2 directory

I want to create a simple snailplot, to understand how it works.
First, I create the necessary files to produce plots afterwards.

``` bash
Sinteractive -c 5 -m 50G -t 01:00:00
# Activate "btk" conda environment

# Father
blobtools create --fasta ../0_FinalFiles/Final_FatherAssembly.fa --busco ./3_BUSCO/BUSCO_Final_FatherAssembly.fa/run_aves_odb10/full_table.tsv --taxid 56313 --taxdump ./2_NCBI TytoAlba_Father_BlobDir
# Mother
blobtools create --fasta ../0_FinalFiles/Final_MotherAssembly.fa --busco ./3_BUSCO/BUSCO_Final_MotherAssembly.fa/run_aves_odb10/full_table.tsv --taxid 56313 --taxdump ./2_NCBI TytoAlba_Mother_BlobDir
```

## Snail plot

From the files in `TytoAlba_*_BlobDir`, I can create a minimal snailplot.

``` bash
# Activate "btk" conda environment

# Father
blobtools view --plot --view snail --format svg TytoAlba_Father_BlobDir
# Mother
blobtools view --plot --view snail --format svg TytoAlba_Mother_BlobDir
```

## UNZIP

⚠️⚠️⚠️ This analysis is not used in the current version of the paper. ⚠️⚠️⚠️

To compare what we have now with what we had previously, I assess the BUSCO score of the new annotated set of genes of the UNZIP and will do the same with the previous version of it (NCBI annotation).

``` bash
mkdir -p 4_UNZIP
# New UNZIP
sbatch 1_BUSCO.slurm /work/FAC/FBM/DEE/jgoudet/barn_owl/hcorval/Chapter_NewReferenceGenome/13_GenomeAnnotation/4_UNZIPAnnotation/UNZIPAnnotation/complete.proteins.faa NewUNZIP
sbatch 1_BUSCO.slurm GCF_018691265.1_T.alba_DEE_v4.0_protein.faa OldUNZIP
```

## Quast

With `Quast`, I should be able to compute the general assembly metrics for a given fasta.

``` bash
# Installation
wget https://github.com/ablab/quast/releases/download/quast_5.3.0/quast-5.3.0.tar.gz
tar -xzf quast-5.3.0.tar.gz
cd quast-5.3.0
```

I can run the Python script using only the fasta files to get the assembly metrics.

``` bash
# Father
python3 quast-5.3.0/quast.py ../0_FinalFiles/Final_FatherAssembly.fa -o 5.1_QuastFather
# Mother
python3 quast-5.3.0/quast.py ../0_FinalFiles/Final_MotherAssembly.fa -o 5.2_QuastMother
```

I also need to compute the general assembly metrics on the non-scaffolded unitigs. I copy-pasted them in my `6_UnscaffoldedUnitigs` directory and used `Quast` to do the rest.

``` bash
mkdir -p 6_UnscaffoldedUnitigs
# Father
cp /nas/FAC/FBM/DEE/jgoudet/barn_owl/D2c/bix/hifiasm_ctg/HiFi/Fa.0518.p_ctg.fa.gz 6_UnscaffoldedUnitigs/
python3 quast-5.3.0/quast.py 6_UnscaffoldedUnitigs/Fa.0518.p_ctg.fa.gz -o 6.1_QuastFather_Unscaffolded
# Mother
cp /nas/FAC/FBM/DEE/jgoudet/barn_owl/D2c/bix/hifiasm_ctg/HiFi/Mo.0518.p_ctg.fa.gz 6_UnscaffoldedUnitigs/
python3 quast-5.3.0/quast.py 6_UnscaffoldedUnitigs/Mo.0518.p_ctg.fa.gz -o 6.2_QuastMother_Unscaffolded
```

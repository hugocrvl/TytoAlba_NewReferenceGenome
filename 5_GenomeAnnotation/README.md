
# Genome annotation

## Table of contents

- [Description](#description)
- [Example](#launching-an-example-of-annotation)
- [Father haplotype](#father-haplotype)
  - [RNA-seq data](#rna-seq-data)
  - [EGAPx configuration](#egapx-configuration)
  - [Annotation](#annotation)
- [UNZIP annotation](#unzip-annotation)
- [Mother annotation](#mother-haplotype)
- [Exploring output](#exploring-output)

## Description

In this directory, I plan to annotate the new genome we have. Three main approaches are possible: _ad initio_, homology-based and RNA-seq. For each approach, a lot of different tools can be used, so I decided to use the publicly available NCBI pipeline for eukaryotic genomes: Eukaryotic Genome Annotation Pipeline External (EGAPx).

## Launching an example of annotation

Now that I've configured my files, I can try launching an example of annotation to see if everything works.

``` bash
# Starting a screen for the pipeline not to close before the end
screen -r EGAPx
# Loading required modules: Nextflow, Java17, Singularity
module load nextflow
module load openjdk
module load singularityce
# Activating the correct virtual environment
source egapx/venv/bin/activate
# Launching the pipeline
python egapx/ui/egapx.py egapx/examples/input_D_farinae_small.yaml -e slurm -c egapx/egapx_config/ -o Example_Out
```

After, I erased the `Example_Out` directory.

## Father haplotype

### RNA-seq data

I need to gather all the RNA-seq reads, copy them into my work, get their tag and paths put them inside a `.txt` file for the pipeline to use it.

``` bash
mkdir -p 2_RNAseq
# 1 individuals, 6 different tissues, keep only paired-end reads
shopt -s extglob
cp /nas/FAC/FBM/DEE/jgoudet/barn_owl/D2c/tcumer/SCRATCH_AXIOM/Analyses_Genome/Analyses_Unzip/4_RNAseq/Reads/!(lib3_*)_paired.fq 2_RNAseq/
# Creating the txt file with all paths
ls 2_RNAseq/ | cut -f1,2 -d_ > 2.1_Path_RNAseq.txt
paste <(ls 2_RNAseq/ | cut -d_ -f1,2) <(realpath 2_RNAseq/*) > 2.1_Path_RNAseq.txt
```

### EGAPx configuration

I need to create directory for Netxflow to write its files inside and specifiy these paths inside the `egapx/egapx_config/slurm.config` file.

``` bash
mkdir -p 3.1_OutNextflow/singularity
mkdir -p 3.1_OutNextflow/tmp
```

I also had to create the `.yaml` config file with `FASTA` path, `taxid`, `RNA-seq` path.

### Annotation

Now that I've seen how it works, I can launch the annotation of our genome the same way as the example.

``` bash
# Starting a screen for the pipeline not to close before the end
screen -S EGAPx_Father
# Loading required modules: Nextflow, Java17, Singularity
module load nextflow
module load openjdk
module load singularityce
# Activating the correct virtual environment
source egapx/venv/bin/activate
# Launching the pipeline
python egapx/ui/egapx.py 1.1_Father_TytoAlba_Configuration.yaml -e slurm -c ./egapx/egapx_config/ -w ./3.1_OutNextflow/workingcloudexecutor -o 3_FatherAnnotation > 3.2_FatherAnnotation.log 2>&1
```

## UNZIP annotation

If we want to quantify the annotation improvment, we need to have the UNZIP assembly annotated the same way as the new reference genome. For this, I'm running the same EGAPx pipeline using the UNZIP assembly but with another `.yaml` and `slurm.config` files

``` bash
# I create a 4_UNZIP directory and work from it to not mix outputs
mkdir -p 4_UNZIPAnnotation
cd 4_UNZIPAnnotation
# I create the singularity and tmp directories for the pipeline to write files inside and specifiy its path in a new slurm.config file
mkdir -p singularity
mkdir -p tmp
mkdir -p egapx_config_UNZIP
cp ../egapx/egapx_config/ egapx_config_UNZIP
# Starting a screen for the pipeline for it not to close before the end
screen -S EGAPx_UNZIP
# Loading required modules: Nextflow, Java17, Singularity
module load nextflow
module load openjdk
module load singularityce
# Activating the correct virtual environment
source ../egapx/venv/bin/activate
# Launching the pipeline
python ../egapx/ui/egapx.py 1_UNZIP_TytoAlba_Configuration.yaml -e slurm -c ./egapx_config_UNZIP/ -w ./workingcloudexecutor -o UNZIPAnnotation > UNZIPAnnotation.log 2>&1
```

## Mother haplotype

After consideration, we'd like to annotate the Mother haplotype as well. I just need to create a specific `.yaml` file with the other reference, paths of the `singularity` and `tmp` directories in the `slurm.config` file as well as the output directory.

``` bash
# Creation of the directories
mkdir -p 5.1_OutNextflow/singularity # Need to change the path in slurm.config
mkdir -p 5.1_OutNextflow/tmp # Need to change the path in slurm.config
# Starting a screen for the pipeline not to close before the end
screen -S EGAPx_Mother
# Loading required modules: Nextflow, Java17, Singularity
module load nextflow
module load openjdk
module load singularityce
# Activating the correct virtual environment
source egapx/venv/bin/activate
# Launching the pipeline
python egapx/ui/egapx.py 1.2_Mother_TytoAlba_Configuration.yaml -e slurm -c ./egapx/egapx_config/ -w ./5.1_OutNextflow/workingcloudexecutor -o 5_MotherAnnotation > 5.2_MotherAnnotation.log 2>&1
```

## Exploring output

The annotation pipeline resulted in a lot of different output files. I want to quantify the amount of genes per windows of 100-kbp and use it for the circos plot. Everything is contained in `3_FatherAnnotation/complete.genomic.gff` and `5_FatherAnnotation/complete.genomic.gff`.

### Comparison with UNZIP

We want to know which genes are new compared to the UNZIP annotation. The way we found was to map with `minimap2` the annotated Coding DNA sequences (CDS) against the UNZIP assembly, and keep CDS that are perfectly aligned, full length. Those that we filter out can be considered as "new" compared to the UNZIP assembly.

However, we don't want to map every single genes to the UNZIP assembly, but only the ones that are also fully mapped, with MAPQ = 60 on their respective haplotype. Therefore, I first map Father CDS to the Father haplotype, and Mother CDS to the Mother haplotype, filter for the best aligned CDS and align only these CDS to the UNZIP.

``` bash
mkdir -p 7.2_AlignedCDS2ASSEMBLY
module load seqkit
# Father CDS to Father assembly
sbatch 7_MappingCDS2ASSEMBLY.slurm ../0_FinalFiles/Final_FatherAssembly.fa 3_FatherAnnotation/complete.cds.fna 7.2_AlignedCDS2ASSEMBLY/FatherCDS_Aligned2FatherAssembly.paf
awk '$12 == 60 && $4 - $3 == $2 {print $0}' 7.2_AlignedCDS2ASSEMBLY/FatherCDS_Aligned2FatherAssembly.paf > 7.2_AlignedCDS2ASSEMBLY/FatherCDS_Aligned2FatherAssembly_Filtered.paf
awk '{print $1}' 7.2_AlignedCDS2ASSEMBLY/FatherCDS_Aligned2FatherAssembly_Filtered.paf > 7.2_AlignedCDS2ASSEMBLY/FatherCDS_Aligned2FatherAssembly_Filtered_NamedOnly.txt
seqkit grep -f 7.2_AlignedCDS2ASSEMBLY/FatherCDS_Aligned2FatherAssembly_Filtered_NamedOnly.txt 3_FatherAnnotation/complete.cds.fna > 3_FatherAnnotation/filtered.cds.fna
# Mother CDS to Mother assembly 
sbatch 7_MappingCDS2ASSEMBLY.slurm ../0_FinalFiles/Final_MotherAssembly.fa 5_MotherAnnotation/complete.cds.fna 7.2_AlignedCDS2ASSEMBLY/MotherCDS_Aligned2MotherAssembly.paf
awk '$12 == 60 && $4 - $3 == $2 {print $0}' 7.2_AlignedCDS2ASSEMBLY/MotherCDS_Aligned2MotherAssembly.paf > 7.2_AlignedCDS2ASSEMBLY/MotherCDS_Aligned2MotherAssembly_Filtered.paf
awk '{print $1}' 7.2_AlignedCDS2ASSEMBLY/MotherCDS_Aligned2MotherAssembly_Filtered.paf > 7.2_AlignedCDS2ASSEMBLY/MotherCDS_Aligned2MotherAssembly_Filtered_NamedOnly.txt
seqkit grep -f 7.2_AlignedCDS2ASSEMBLY/MotherCDS_Aligned2MotherAssembly_Filtered_NamedOnly.txt 5_MotherAnnotation/complete.cds.fna > 5_MotherAnnotation/filtered.cds.fna
```

Now, I have the CDS from both haplotypes (`*_*Annotation/filtered.cds.fna`) that I should map to the UNZIP assembly, then filter for full length and MAPQ = 60.

``` bash
# Father
sbatch 7_MappingCDS2ASSEMBLY.slurm ../../Chapter_FitnessPhenotypicEffects/0_ReferenceGenome/Tyto_reference_Jan2020.fasta 3_FatherAnnotation/filtered.cds.fna 7.2_AlignedCDS2ASSEMBLY/FilteredFatherCDS_Aligned2UNZIP.paf
awk '{print $1, $2, $3, $4, $6, $12}' 7.2_AlignedCDS2ASSEMBLY/FilteredFatherCDS_Aligned2UNZIP.paf > 7.2_AlignedCDS2ASSEMBLY/FilteredFatherCDS_Aligned2UNZIP_Subset.txt
# Mother
sbatch 7_MappingCDS2ASSEMBLY.slurm ../../Chapter_FitnessPhenotypicEffects/0_ReferenceGenome/Tyto_reference_Jan2020.fasta 5_MotherAnnotation/filtered.cds.fna 7.2_AlignedCDS2ASSEMBLY/FilteredMotherCDS_Aligned2UNZIP.paf
awk '{print $1, $2, $3, $4, $6, $12}' 7.2_AlignedCDS2ASSEMBLY/FilteredMotherCDS_Aligned2UNZIP.paf > 7.2_AlignedCDS2ASSEMBLY/FilteredMotherCDS_Aligned2UNZIP_Subset.txt
```

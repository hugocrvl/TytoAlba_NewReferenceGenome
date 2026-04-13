# Tyto alba – New Reference Genome Pipeline

This repository contains all scripts, workflows, and documentation associated with the generation and analysis of the **new chromosome‑level, haplotype‑resolved reference genome of the Western barn owl (*Tyto alba*)**.

## Repository Structure

Each of these directory has its own README file reporting the commands line to run the steps described in this directory. 

#### 1_AlleleSorting
Tools for separating parental haplotypes.

#### 2_assembly
Scripts and workflows for assembling raw sequencing data.

#### 2_microChromsomes
Identification of micro chromsomes. 

#### 3_QC
Assembly quality metrics and BUSCO analysis.

#### 4_PhasingEvaluation
Assessment of haplotype phasing consistency.

#### 4_Synteny
Comparative genomics and synteny visualization of the two haplotypes.

#### 5_GenomeAnnotation 
Automated and manual annotation steps.

#### 5_RepeatsRed
RepeatModeler/RepeatMasker workflows.

#### 5_RepeatsRepeatMasker
RepeatModeler/RepeatMasker workflows.

#### 6_MappabilityMask
Generation of genomic masks for alignment mappability.

#### 6_RecombinationMapLiftOver
Tools to transfer recombination maps onto the new reference.

#### 7_Circos
Scripts for Circos genomic visualization.
In this directory is also described the repeat detection using RepeatModeler and RepeatMasker.

#### 7_ChromosomeCategories 
Structural chromosome classification tools.

####  8_protein_based_synteny
Comparative genomics and synteny with other bird species


## Citation
If you use this repository, please cite the publication (DOI: 10.64898/2026.03.20.713190).

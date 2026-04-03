
# Circos plot construction

## Table of contents

- [Description](#description)
- [Karyotypes and ideograms](#karyotypes-and-ideograms)
- [Telomeric sequences](#telomeric-sequences)
- [Repeats quantification](#repeats-quantification)
- [Repeats distribution](#repeats-distribution)
- [Alignment and ribbon links between chromosomes](#ribbon-links-between-chromosomes)
- [Gene density across the genome](#gene-density-across-the-genome)
- [Circos creation](#circos-creation)

## Description

In this directory, I build input files and what is necessary for the construction of a circos plot for the new reference genome paper. We plan to include:

- Paternal and maternal chromosomes,
- telomeric and centromeric locations if available.
- a central alignment between the two haplotypes,
- repeats distribution,
- a mappability profil,
- gene density,

## Circos installation

Circos is a set of tools that can be installed by doing:

``` bash
module load miniforge3
conda_init
conda create -n circos -c conda-forge perl
conda activate circos
conda install -c bioconda circos
circos --help
conda deactivate
```

To try making a mock figure:

``` bash
# Activate circos conda environment
circos -conf ~/WORK/CONDA/conda.environments/circos/example/etc/circos.conf
```

The `circos.conf` file is the mandatory and backbone file of `circos`.

## Karyotypes and ideograms

To create the karyotypes for each haplotype, run:

``` bash
# Creation of an initial file for the Father chromosomes
bash ./1.1_KaryotypeConstruction.sh ../0_FinalFiles/Final_FatherAssembly.fa chrfather > ./TextFiles/Father_Karyotype.txt

# Adding the Mother chromosomes
bash ./1.1_KaryotypeConstruction.sh ../0_FinalFiles/Final_MotherAssembly.fa chrmother > ./TextFiles/Mother_Karyotype.txt
```

I had to define the colours I'll use for the different chromosomes beforehand, in `etc/colors.conf`. I added `chrfather = 21, 96, 121` and `chrmother = 136, 90, 90`.

To specifiy the order of the chromosomes, I added the `chromosomes_order` parameter inside the `./ConfigurationFiles/circos.conf` file.
Another important parameter is `chromosomes_units` to define the size of unit, useful to add space between chromosomes.

To display the chromosome labels, I added several parameters in the `./ConfigurationFiles/circos.conf` file, see `<spacing>` code chunck.

To change the png name, offset or orientation angle, look at the files : `./ConfigurationFiles/image.conf` > `./ConfigurationFiles/image.generic.conf`.

## Telomeric sequences

For every chromosome for which the telomeres is known, I add a corresponding circle on the ideogram/chromosomes.

To do so, I manually created a `.TextFiles/tmp` made of the telomeric positions from `/nas/FAC/FBM/DEE/jgoudet/barn_owl/D2c/bix/202512_final/all_telo.xlsx`. I merged the windows that were next to each other, to have maximum two entries (= 2 telomeres) per chromosome. Then, I renamed the chromosomes according to the last version of the genome and added the value `1` with the `1.2_ConstructionTelomericFile.sh` script.

I also added `telomeres = 255, 177, 0` in the same file as before, `etc/colors.conf`.

``` bash
bash 1.2_ConstructionTelomericFile.sh TextFiles/tmp.txt ../0_FinalFiles/1_SortingRenaming/2.1_CorrespondanceFather.txt ../0_FinalFiles/1_SortingRenaming/2.2_CorrespondanceMother.txt 1 > ./TextFiles/Telomeric_Positions.txt
rm ./TextFiles/tmp.txt
```

To display the telomeric positions, I added several parameters in the `./ConfigurationFiles/circos.conf` file, see `<plots>` code chunck.

## Repeats quantification

I quantify the amount of repeats per assembly. This will be displayed as an outer layer of the circos plot, but also regions of high repeat density will be masked in the assembly alignment. I use `RepeatModeler` to characterise the repeats, and `RepeatMasker` to detect and mask them.

First, I need to create a database from my `fasta` files, one per assembly.

``` bash
module load repeatmodeler

mkdir -p 2_Repeats
cd 2_Repeats

# Father
mkdir -p Father
cd Father
BuildDatabase -name TytoDB ../../../0_FinalFiles/Final_FatherAssembly.fa

# Mother
mkdir -p Mother # in 2_Repeats
cd Mother
BuildDatabase -name TytoDB ../../../0_FinalFiles/Final_MotherAssembly.fa
```

Then, I run `RepeatModeler`. It will detect repeats.

``` bash
# Father
sbatch 2.1_RunRepeatModeler.slurm /work/FAC/FBM/DEE/jgoudet/barn_owl/hcorval/Chapter_NewReferenceGenome/0_FinalFiles/Final_FatherAssembly.fa ./2_Repeats/Father/

# Mother
sbatch 2.1_RunRepeatModeler.slurm /work/FAC/FBM/DEE/jgoudet/barn_owl/hcorval/Chapter_NewReferenceGenome/0_FinalFiles/Final_MotherAssembly.fa ./2_Repeats/Mother/
```

And `RepeatMasker`.

I need to wait for the Father pipeline to finish before running the Mother one, as it uses the same cache file: `/users/hcorval/.RepeatMaskerCache//general.working/is.lib`.

If I want to change the run of `RepeatModeler`, I can simply change the `RM_*` directory of the input file, as well as the output directory.

``` bash
# Father
mkdir -p 2_Repeats/Father/Out_RepeatMasker
sbatch 2.2_RunRepeatMasker.slurm /work/FAC/FBM/DEE/jgoudet/barn_owl/hcorval/Chapter_NewReferenceGenome/0_FinalFiles/Final_FatherAssembly.fa 2_Repeats/Father/RM_2419703.FriFeb202232422026/consensi.fa.classified 2_Repeats/Father/Out_RepeatMasker

# Mother
mkdir -p 2_Repeats/Mother/Out_RepeatMasker
sbatch 2.2_RunRepeatMasker.slurm /work/FAC/FBM/DEE/jgoudet/barn_owl/hcorval/Chapter_NewReferenceGenome/0_FinalFiles/Final_MotherAssembly.fa 2_Repeats/Mother/RM_2993976.SatFeb210036482026/consensi.fa.classified 2_Repeats/Mother/Out_RepeatMasker
```

## Repeats distribution

I quantify the percentage of repeats per windows of 100 kbp.

``` bash
# Father
bash ./2.3_QuantifyRepeatDistribution.sh ./2_Repeats/Father/Out_RepeatMasker/Final_FatherAssembly.fa.out > ./TextFiles/RepeatDistribution_Father.txt
# Mother
bash ./2.3_QuantifyRepeatDistribution.sh ./2_Repeats/Mother/Out_RepeatMasker/Final_MotherAssembly.fa.out > ./TextFiles/RepeatDistribution_Mother.txt
```

To display the repeats distribution, I added several parameters in the `./ConfigurationFiles/circos.conf` file, see `<plots>` code chunck.

## Ribbon links between chromosomes

I need to map the Father assembly to the Mother assembly to obtain the whole genome alignment. I do it with `minimap2`, and set the output in the `paf` format.

``` bash
mkdir -p 3.2_WholeGenomeAlignments

# Genome with masked repeats, RepeatModeler and RepeatMasker pipeline
sbatch 3_WholeGenomeAlignment.slurm ./2_Repeats/Father/Out_RepeatMasker/Final_FatherAssembly.fa.masked ./2_Repeats/Mother/Out_RepeatMasker/Final_MotherAssembly.fa.masked 3.2_WholeGenomeAlignments/WholeGenome_MaskedRepeats.paf
```

Then, I need to split the `paf` alignment into two `txt` files: one having only the `+` orientation, the other with the `-` orientation. I keep alignement of 1000 kbp minimum, otherwise everything map everywhere.

``` bash
# Similar orientation
bash 3.3_LinksFilesConstruction.sh ./3.2_WholeGenomeAlignments/WholeGenome_MaskedRepeats.paf + 10000 ./TextFiles/Alignment_SimilarOrientation.txt
# Opposite orientation
bash 3.3_LinksFilesConstruction.sh ./3.2_WholeGenomeAlignments/WholeGenome_MaskedRepeats.paf - 10000 ./TextFiles/Alignment_OppositeOrientation.txt
```

To add the links, I add `<links>` code chuncks, one per alignment orientation, in the `./ConfigurationFiles/circos.conf` file.

I also added additional colours in `etc/colors.conf`. I added `linkssimilar = 128, 143, 135` and `linksopposite = 67, 76, 71`.

## Gene density across the genome

I've annotated the genome using the EGAPx pipeline. Analysis can be found in `../13_GenomeAnnotation` directory. Now, I want to quantify the amount of genes per haplotype and I can use a customed script for this.

``` bash
# Father
bash 5_GeneDensity/1_GD_Quantification.sh ../13_GenomeAnnotation/6_ProcessingOutput/1.1_GenomicAnnotation_Fa_GenesOnly.txt 100000 > ./TextFiles/GD_Father.txt
# Mother
bash 5_GeneDensity/1_GD_Quantification.sh ../13_GenomeAnnotation/6_ProcessingOutput/1.2_GenomicAnnotation_Mo_GenesOnly.txt 100000 > ./TextFiles/GD_Mother.txt
```

## Circos creation

To obtain a figure with my `.conf` file, run:

``` bash
# Activate circos conda environment
circos -conf ./ConfigurationFiles/circos.conf
```

## Split HiFi reads in categories based on sequence features

```console
./splitFqHifiReads
```

## Futher split HiFi reads based on expected allelic origin

```console
./gatherMergedFaMo5
```

## Version of hifiasm used for HiFi reads assembly

```console
hifiasm  --version
0.19.8-r603
```

## Assembly of Father and Mother alleles

```console
hifiasm -o /data/chris/tyto/20240518Fa/asm -t 72 --primary -D 10.0 -N 1000 --max-kocc 20000 --n-hap 1 --ul SOL000x_SUP_20k_cutdup.fq.gz,LR_40k.fq.gz /data/chris/tyto/20240518Fa/{other,musat,telomers}.fq.gz /data/chris/tyto/20240518Fa/Unslct.{other,musat,telomers}.fq.gz 2>&1 | tee 20240518Fa_log.txt
hifiasm -o /data/chris/tyto/20240518Mo/asm -t 72 --primary -D 10.0 -N 1000 --max-kocc 20000 --n-hap 1 --ul SOL000x_SUP_20k_cutdup.fq.gz,LR_40k.fq.gz /data/chris/tyto/20240518Mo/{other,musat,telomers}.fq.gz /data/chris/tyto/20240518Fa/Unslct.{other,musat,telomers}.fq.gz 2>&1 | tee 20240518Mo_log.txt
```

## Mapping of some mRNA of interest and all HiFi reads to the assembled contigs

```console
pushd /data/chris/tyto/20240518Fa/
minimap2 -x splice:hq asm.p_ctg.fa ~/tyto/SeqGeneForChristian.fasta >Fa_mRNA.paf
cd ../20240518Mo/
minimap2 -x splice:hq asm.p_ctg.fa ~/tyto/SeqGeneForChristian.fasta >Mo_mRNA.paf

./doMapMoFa 
```

## Compute coverage by the HiFi reads on the assembled contigs

```console
sort -m -k 6,6 -k 8,8n -S 10% /data/chris/tyto/20240518Fa/m64???_*.paf >m64_all_Fa.paf &
sort -m -k 6,6 -k 8,8n -S 10% /data/chris/tyto/20240518Mo/m64???_*.paf >m64_all_Mo.paf &

./genCoverage m64_all_Fa.paf >m64_all_Fa_cov.txt &
./genCoverage m64_all_Mo.paf >m64_all_Mo_cov.txt &

./analyzeCov >zero_cov_regions.txt
```

## Determine zero coverage regions to be further investigated

```console
./grabCandidates zero_cov_regions.txt >zero_cov_regions_tags.txt
```


## Version of hifiasm used for Oxford Nanopore reads assembly

```console
hifiasm --version
0.25.0-r726
```

## Split ON reads in categories based on sequence features and expected allelic origin

```console
./prepareONreads
```

## Assembly of Father and Mother alleles using ON reads only

```console
hifiasm --ont -o hifiasm.Fa/asm -t 40 --primary --n-hap 1 -l 0 --min-hist-cnt 3 ON?.{Father,Unslct}.{other,musat,telo}.fq.gz 2>&1 | tee hifiasm_Fa_log.txt
hifiasm --ont -o hifiasm.Mo/asm -t 40 --primary --n-hap 1 -l 0 --min-hist-cnt 3 ON?.{Mother,Unslct}.{other,musat,telo}.fq.gz 2>&1 | tee hifiasm_Mo_log.txt

for p in Fa Mo; do grep ^S hifiasm.$p/asm.p_ctg.gfa | perl -ane 'printf ">%s\n%s\n", $F[1], $F[2]' >hifiasm.$p/asm.p_ctg.fa; done

sed 's/^>/>ON_Fa_/' hifiasm.Fa/asm.p_ctg.fa >hifiasm_Fa.fa
sed 's/^>/>ON_Mo_/' hifiasm.Mo/asm.p_ctg.fa >hifiasm_Mo.fa
```

## Mapping of assembled HiFi contigs to the assembled ON contigs

```console
for p1 in Fa Mo; do for p2 in Fa Mo; do echo "$p1 vs $p2"; minimap2 -t 40 -x asm5 -I 12G --secondary=no $p1.0518.p_ctg.fa.gz hifiasm_$p2.fa |perl -ane 'print if ($F[3]-$F[2])>=40000'>$p1.0518.vs.$p2.paf; done; done
```

## Prepare alignments targeting the zero coverage regions and Regions Of Interest between assembled contigs

```console
./runLastz

./doMafft ROI_regions_tags.txt

./move2LG 
```

## Iterative process to create each polished chromosome based on the Linkage Groups

```console
./getAnchors aln/LG01 >LG01_anchors.txt
cp LG01_anchors.txt LG01_build.txt 
vi LG01_build.txt 
./buildIt LG01_build.txt >LG01.fa
head -2 LG01.fa >LG01_Fa.fa 
tail -2 LG01.fa >LG01_Mo.fa 
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG01_Fa.fa LG01_Mo.fa >/tmp/m.paf
minidot /tmp/m.paf >/tmp/m.eps
epstopdf /tmp/m.eps 

./getAnchors aln/LG02 >LG02_anchors.txt
cp LG02_anchors.txt LG02_build.txt 
cat LG01_build.txt >>LG02_build.txt 
vi LG02_build.txt 
./buildIt LG02_build.txt >LG02.fa
head -2 LG02.fa >LG02_Fa.fa 
tail -2 LG02.fa >LG02_Mo.fa 
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG02_Fa.fa LG02_Mo.fa >/tmp/m.paf
minidot /tmp/m.paf >/tmp/m.eps
epstopdf /tmp/m.eps 

./getAnchors aln/LG03 >LG03_anchors.txt
cp LG03_anchors.txt LG03_build.txt 
cat LG02_build.txt >>LG03_build.txt 
vi LG03_build.txt 
./buildIt LG03_build.txt >LG03.fa
head -2 LG03.fa >LG03_Fa.fa 
tail -2 LG03.fa >LG03_Mo.fa 
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG03_Fa.fa LG03_Mo.fa >/tmp/m.paf
minidot /tmp/m.paf >/tmp/m.eps
epstopdf /tmp/m.eps 

./getAnchors aln/LG04 >LG04_anchors.txt
less LG04_anchors.txt 
cp LG04_anchors.txt LG04_build.txt 
cat LG03_build.txt >>LG04_build.txt 
vi LG04_build.txt 
./buildIt LG04_build.txt >LG04.fa
head -2 LG04.fa >LG04_Fa.fa 
tail -2 LG04.fa >LG04_Mo.fa 
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG04_Fa.fa LG04_Mo.fa >/tmp/m.paf
minidot /tmp/m.paf >/tmp/m.eps
epstopdf /tmp/m.eps 
cat LG??_Fa.fa >LG_Fa.fa
cat LG??_Mo.fa >LG_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG_Fa.fa LG_Mo.fa >/tmp/m.paf
minidot /tmp/m.paf >/tmp/m.eps
epstopdf /tmp/m.eps 

./getAnchors aln/LG05 >LG05_anchors.txt
less LG05_anchors.txt 
perl -ane 'print if $F[10] eq "Mo_000019l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 13,13n|less
head -1 /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt
perl -ane 'print if $F[19] eq "Mo_000019l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 22,22n|less
vi doMafft 
./doMafft ROI_regions_tags.txt
rm aln/Fa_000005l_12252186_12287599_all.*
vi ROI_regions_tags.txt 
./doMafft ROI_regions_tags.txt
rm aln/Mo_000019l_31567601_31606358_*
./move2LG 
./getAnchors aln/LG05 >LG05_anchors.txt
less LG05_anchors.txt 
cp LG05_anchors.txt LG05_build.txt 
cat LG05_build.txt >>LG05_build.txt 
cat LG04_build.txt >>LG05_build.txt 
vi LG05_build.txt 
./buildIt LG05_build.txt >LG05.fa
head -2 LG05.fa >LG05_Fa.fa 
tail -2 LG05.fa >LG05_Mo.fa 
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG05_Fa.fa LG05_Mo.fa >/tmp/m.paf
minidot /tmp/m.paf >/tmp/m.eps
epstopdf /tmp/m.eps 

./getAnchors aln/LG06 >LG06_anchors.txt
less LG06_anchors.txt 
perl -ane 'print if $F[19] eq "Mo_000005l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 22,22n|less
./doMafft ROI_regions_tags.txt
./move2LG 
./getAnchors aln/LG06 >LG06_anchors.txt
less LG06_anchors.txt 
cp LG06_anchors.txt LG06_build.txt 
cat LG05_build.txt >>LG06_build.txt 
vi LG06_build.txt 
grep Fa_000132l /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt|head|less
vi LG06_build.txt 
./buildIt LG06_build.txt >LG06.fa
head -2 LG06.fa >LG06_Fa.fa 
tail -2 LG06.fa >LG06_Mo.fa 
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG06_Fa.fa LG06_Mo.fa >/tmp/m.paf
minidot /tmp/m.paf >/tmp/m.eps
epstopdf /tmp/m.eps 
cat LG??_Fa.fa >LG_Fa.fa
cat LG??_Mo.fa >LG_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG_Fa.fa LG_Mo.fa >/tmp/m.paf
minidot /tmp/m.paf >/tmp/m.eps
epstopdf /tmp/m.eps 

cp LG06_build.txt LG07_build.txt 
vi LG07_build.txt 
grep Mo_000049l /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt|head -2 >>LG07_build.txt 
grep Fa_000015l /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt|head -2 >>LG07_build.txt 
vi LG07_build.txt 
./buildIt LG07_build.txt >LG07.fa
head -2 LG07.fa >LG07_Fa.fa 
tail -2 LG07.fa >LG07_Mo.fa 
cat LG??_Fa.fa >LG_Fa.fa
cat LG??_Mo.fa >LG_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG_Fa.fa LG_Mo.fa >/tmp/m.paf
minidot /tmp/m.paf >/tmp/m.eps
epstopdf /tmp/m.eps 

cp LG07_build.txt LG08_build.txt 
grep Mo_000017l /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt|head -2 >>LG08_build.txt 
vi LG08_build.txt 
./buildIt LG08_build.txt >LG08.fa
head -2 LG08.fa >LG08_Fa.fa 
tail -2 LG08.fa >LG08_Mo.fa 
cat LG??_Fa.fa >LG_Fa.fa
cat LG??_Mo.fa >LG_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG_Fa.fa LG_Mo.fa >/tmp/m.paf
minidot /tmp/m.paf >/tmp/m.eps
epstopdf /tmp/m.eps 

./getAnchors aln/LG09 >LG09_anchors.txt
less LG09_anchors.txt 
perl -ane 'print if $F[19] eq "Mo_000075l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 22,22n|less
./doMafft ROI_regions_tags.txt
./move2LG 
./getAnchors aln/LG09 >LG09_anchors.txt
less LG09_anchors.txt 
cp LG09_anchors.txt LG09_build.txt 
cat LG08_build.txt >>LG09_build.txt 
vi LG09_build.txt 
./buildIt LG09_build.txt >LG09.fa
head -2 LG09.fa >LG09_Fa.fa 
tail -2 LG09.fa >LG09_Mo.fa 
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG09_Fa.fa LG09_Mo.fa >/tmp/m.paf
minidot /tmp/m.paf >/tmp/m.eps
epstopdf /tmp/m.eps 

./getAnchors aln/LG10 >LG10_anchors.txt
less LG10_anchors.txt 
perl -ane 'print if $F[10] eq "Fa_000105l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 13,13n|less
./doMafft ROI_regions_tags.txt
./move2LG 
./getAnchors aln/LG10 >LG10_anchors.txt
less LG10_anchors.txt 
ls aln/LG10/
less aln/LG10/Fa_000105l_22927074_22938956_log.txt 
cp LG09_build.txt LG10_build.txt 
vi LG10_build.txt 
cat LG10_anchors.txt >>LG10_build.txt 
vi LG10_build.txt 
./buildIt LG10_build.txt >LG10.fa
head -2 LG10.fa >LG10_Fa.fa 
tail -2 LG10.fa >LG10_Mo.fa 
cat LG??_Fa.fa >LG_Fa.fa
cat LG??_Mo.fa >LG_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG_Fa.fa LG_Mo.fa >/tmp/m.paf
minidot /tmp/m.paf >/tmp/m.eps
epstopdf /tmp/m.eps 
vi LG10_build.txt 
./buildIt LG10_build.txt >LG10.fa
head -2 LG10.fa >LG10_Fa.fa 
tail -2 LG10.fa >LG10_Mo.fa 
cat LG??_Fa.fa >LG_Fa.fa
cat LG??_Mo.fa >LG_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG_Fa.fa LG_Mo.fa >/tmp/m.paf
minidot /tmp/m.paf >/tmp/m.eps
epstopdf /tmp/m.eps 

./getAnchors aln/LG11 >LG11_anchors.txt
less LG11_anchors.txt 
cp LG11_anchors.txt LG11_build.txt 
cat LG10_build.txt >>LG11_build.txt 
vi LG11_build.txt 
./buildIt LG11_build.txt >LG11.fa
head -2 LG11.fa >LG11_Fa.fa 
tail -2 LG11.fa >LG11_Mo.fa 
cat LG??_Fa.fa >LG_Fa.fa
cat LG??_Mo.fa >LG_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG_Fa.fa LG_Mo.fa >/tmp/m.paf
minidot /tmp/m.paf >/tmp/m.eps
epstopdf /tmp/m.eps 

./getAnchors aln/LG12 >LG12_anchors.txt
less LG12_anchors.txt 
cp LG12_anchors.txt LG12_build.txt 
cat LG11_build.txt >>LG12_build.txt 
vi LG12_build.txt 
./buildIt LG12_build.txt >LG12.fa
cat LG??_Fa.fa >LG_Fa.fa
cat LG??_Mo.fa >LG_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG_Fa.fa LG_Mo.fa >/tmp/m.paf
minidot /tmp/m.paf >/tmp/m.eps
epstopdf /tmp/m.eps 
head -2 LG12.fa >LG12_Fa.fa 
tail -2 LG12.fa >LG12_Mo.fa 
cat LG??_Mo.fa >LG_Mo.fa
cat LG??_Fa.fa >LG_Fa.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG_Fa.fa LG_Mo.fa >/tmp/m.paf
minidot /tmp/m.paf >/tmp/m.eps
epstopdf /tmp/m.eps 

cp LG07_build.txt LG13_build.txt 
vi LG13_build.txt 
./buildIt LG13_build.txt >LG13.fa
head -2 LG13.fa >LG13_Fa.fa 
tail -2 LG13.fa >LG13_Mo.fa 
cat LG??_Fa.fa >LG_Fa.fa
cat LG??_Mo.fa >LG_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG_Fa.fa LG_Mo.fa >/tmp/m.paf
minidot /tmp/m.paf >/tmp/m.eps
epstopdf /tmp/m.eps 
cat LG0?_Fa.fa >LG0_Fa.fa
cat LG0?_Mo.fa >LG0_Mo.fa
cat LG1?_Fa.fa >LG1_Fa.fa
cat LG1?_Mo.fa >LG1_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG0_Fa.fa LG0_Mo.fa >/tmp/m0.paf
minidot /tmp/m0.paf >/tmp/m0.eps
epstopdf /tmp/m0.eps 
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG1_Fa.fa LG1_Mo.fa >/tmp/m1.paf
minidot /tmp/m1.paf >/tmp/m1.eps
epstopdf /tmp/m1.eps 

cp LG13_build.txt LG14_build.txt 
vi LG14_build.txt 
./buildIt LG14_build.txt >LG14.fa
head -2 LG14.fa >LG14_Fa.fa 
tail -2 LG14.fa >LG14_Mo.fa 
cat LG1?_Fa.fa >LG1_Fa.fa
cat LG1?_Mo.fa >LG1_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG1_Fa.fa LG1_Mo.fa >/tmp/m1.paf
minidot /tmp/m1.paf >/tmp/m1.eps
epstopdf /tmp/m1.eps 

./getAnchors aln/LG15 >LG15_anchors.txt
less LG15_anchors.txt 
cp LG15_anchors.txt LG15_build.txt 
cat LG14_build.txt >>LG15_build.txt 
vi LG15_build.txt 
./buildIt LG15_build.txt >LG15.fa
cat LG1?_Fa.fa >LG1_Fa.fa
cat LG1?_Mo.fa >LG1_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG1_Fa.fa LG1_Mo.fa >/tmp/m1.paf
minidot /tmp/m1.paf >/tmp/m1.eps
epstopdf /tmp/m1.eps 
head -2 LG15.fa >LG15_Fa.fa 
tail -2 LG15.fa >LG15_Mo.fa 
cat LG1?_Fa.fa >LG1_Fa.fa
cat LG1?_Mo.fa >LG1_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG1_Fa.fa LG1_Mo.fa >/tmp/m1.paf
minidot /tmp/m1.paf >/tmp/m1.eps
epstopdf /tmp/m1.eps 

./getAnchors aln/LG16 >LG16_anchors.txt
perl -ane 'print if $F[19] eq "Mo_000001l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 22,22n|less
./doMafft ROI_regions_tags.txt
./move2LG 
./getAnchors aln/LG16 >LG16_anchors.txt
less LG16_anchors.txt 
cp LG16_anchors.txt LG16_build.txt 
cat LG15_build.txt >>LG16_build.txt 
vi LG16_build.txt 
./buildIt LG16_build.txt >LG16.fa
head -2 LG65.fa >LG16_Fa.fa 
head -2 LG16.fa >LG16_Fa.fa 
tail -2 LG16.fa >LG16_Mo.fa 
cat LG1?_Fa.fa >LG1_Fa.fa
cat LG1?_Mo.fa >LG1_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG1_Fa.fa LG1_Mo.fa >/tmp/m1.paf
minidot /tmp/m1.paf >/tmp/m1.eps
epstopdf /tmp/m1.eps 
vi /scratch/tyto/plot_cM_pos_11_20.R 
vi LG16_build.txt 
./buildIt LG16_build.txt >LG16.fa
head -2 LG16.fa >LG16_Fa.fa 
tail -2 LG16.fa >LG16_Mo.fa 
cat LG1?_Fa.fa >LG1_Fa.fa
cat LG1?_Mo.fa >LG1_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG1_Fa.fa LG1_Mo.fa >/tmp/m1.paf
minidot /tmp/m1.paf >/tmp/m1.eps
epstopdf /tmp/m1.eps 

./getAnchors aln/LG17 >LG17_anchors.txt
less LG17_anchors.txt 
perl -ane 'print if $F[10] eq "Fa_000057l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 13,13n|less
./doMafft ROI_regions_tags.txt
./move2LG 
./getAnchors aln/LG17 >LG17_anchors.txt
cp LG17_anchors.txt LG17_build.txt 
cat LG16_build.txt >>LG17_build.txt 
vi LG17_build.txt 
./buildIt LG17_build.txt >LG17.fa
head -2 LG17.fa >LG17_Fa.fa 
tail -2 LG17.fa >LG17_Mo.fa 
cat LG1?_Fa.fa >LG1_Fa.fa
cat LG1?_Mo.fa >LG1_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG1_Fa.fa LG1_Mo.fa >/tmp/m1.paf
minidot /tmp/m1.paf >/tmp/m1.eps
epstopdf /tmp/m1.eps 
vi /scratch/tyto/plot_cM_pos_11_20.R 
vi LG17_build.txt 
./buildIt LG17_build.txt >LG17.fa
head -2 LG17.fa >LG17_Fa.fa 
tail -2 LG17.fa >LG17_Mo.fa 
cat LG1?_Fa.fa >LG1_Fa.fa
cat LG1?_Mo.fa >LG1_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG1_Fa.fa LG1_Mo.fa >/tmp/m1.paf
minidot /tmp/m1.paf >/tmp/m1.eps
epstopdf /tmp/m1.eps 

cp LG17_build.txt LG18_build.txt 
vi LG18_build.txt 
perl -ane 'print if $F[19] eq "Mo_000003l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 22,22n|less
ls -R|grep Fa_000091
perl -ane 'print if $F[19] eq "Mo_000003l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 22,22n|less
vi LG18_build.txt 
./buildIt LG18_build.txt >LG18.fa
head -2 LG18.fa >LG18_Fa.fa 
tail -2 LG18.fa >LG18_Mo.fa 
cat LG1?_Fa.fa >LG1_Fa.fa
cat LG1?_Mo.fa >LG1_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG1_Fa.fa LG1_Mo.fa >/tmp/m1.paf
minidot /tmp/m1.paf >/tmp/m1.eps
epstopdf /tmp/m1.eps 

cp LG18_build.txt LG19_build.txt 
vi LG19_build.txt 
./buildIt LG19_build.txt >LG19.fa
head -2 LG19.fa >LG19_Fa.fa 
tail -2 LG19.fa >LG19_Mo.fa 
cat LG1?_Fa.fa >LG1_Fa.fa
cat LG1?_Mo.fa >LG1_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG1_Fa.fa LG1_Mo.fa >/tmp/m1.paf
minidot /tmp/m1.paf >/tmp/m1.eps
epstopdf /tmp/m1.eps 

cp LG19_build.txt LG20_build.txt 
vi LG20_build.txt 
./buildIt LG20_build.txt >LG20.fa
head -2 LG20.fa >LG20_Fa.fa 
tail -2 LG20.fa >LG20_Mo.fa 
cat LG2?_Fa.fa >LG2_Fa.fa
cat LG2?_Mo.fa >LG2_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG2_Fa.fa LG2_Mo.fa >/tmp/m2.paf
minidot /tmp/m1.paf >/tmp/m2.eps
minidot /tmp/m2.paf >/tmp/m2.eps
epstopdf /tmp/m2.eps 

cp LG20_build.txt LG21_build.txt 
vi LG21_build.txt 
perl -ane 'print if $F[10] eq "Fa_000053l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 13,13n|less
vi LG21_build.txt 
./buildIt LG21_build.txt >LG21.fa
head -2 LG21.fa >LG21_Fa.fa 
tail -2 LG21.fa >LG21_Mo.fa 
cat LG2?_Fa.fa >LG2_Fa.fa
cat LG2?_Mo.fa >LG2_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG2_Fa.fa LG2_Mo.fa >/tmp/m2.paf
minidot /tmp/m2.paf >/tmp/m2.eps
epstopdf /tmp/m2.eps 

cp LG21_build.txt LG22_build.txt 
vi LG22_build.txt 
./buildIt LG22_build.txt >LG22.fa
head -2 LG22.fa >LG22_Fa.fa 
tail -2 LG22.fa >LG22_Mo.fa 
cat LG2?_Fa.fa >LG2_Fa.fa
cat LG2?_Mo.fa >LG2_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG2_Fa.fa LG2_Mo.fa >/tmp/m2.paf
minidot /tmp/m2.paf >/tmp/m2.eps
epstopdf /tmp/m2.eps 

./getAnchors aln/LG23 >LG23_anchors.txt
perl -ane 'print if $F[10] eq "Fa_000167l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 13,13n|less
./doMafft ROI_regions_tags.txt
./getAnchors aln/LG23 >LG23_anchors.txt
./move2LG 
./getAnchors aln/LG23 >LG23_anchors.txt
less LG23_anchors.txt 
cp LG22_build.txt LG23_build.txt 
vi LG23_build.txt 
./buildIt LG23_build.txt >LG23.fa
head -2 LG23.fa >LG23_Fa.fa 
tail -2 LG23.fa >LG23_Mo.fa 
cat LG2?_Fa.fa >LG2_Fa.fa
cat LG2?_Mo.fa >LG2_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG2_Fa.fa LG2_Mo.fa >/tmp/m2.paf
minidot /tmp/m2.paf >/tmp/m2.eps
epstopdf /tmp/m2.eps 

./getAnchors aln/LG24 >LG24_anchors.txt
less LG24_anchors.txt 
perl -ane 'print if $F[19] eq "Mo_000016l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 22,22n|less
./doMafft ROI_regions_tags.txt
./move2LG 
./getAnchors aln/LG24 >LG24_anchors.txt
less LG24_anchors.txt 
cp LG24_anchors.txt LG24_build.txt 
cat LG23_build.txt >>LG24_build.txt 
vi LG24_build.txt 
./buildIt LG24_build.txt >LG24.fa
head -2 LG24.fa >LG24_Fa.fa 
tail -2 LG24.fa >LG24_Mo.fa 
cat LG2?_Fa.fa >LG2_Fa.fa
cat LG2?_Mo.fa >LG2_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG2_Fa.fa LG2_Mo.fa >/tmp/m2.paf
minidot /tmp/m2.paf >/tmp/m2.eps
epstopdf /tmp/m2.eps 

cp LG24_build.txt LG25_build.txt 
vi LG25_build.txt 
./buildIt LG25_build.txt >LG25.fa
head -2 LG25.fa >LG25_Fa.fa 
tail -2 LG25.fa >LG25_Mo.fa 
cat LG2?_Fa.fa >LG2_Fa.fa
cat LG2?_Mo.fa >LG2_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG2_Fa.fa LG2_Mo.fa >/tmp/m2.paf
minidot /tmp/m2.paf >/tmp/m2.eps
epstopdf /tmp/m2.eps 

./getAnchors aln/LG26 >LG26_anchors.txt
perl -ane 'print if $F[10] eq "Fa_000065l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 13,13n|less
perl -ane 'print if $F[19] eq "Mo_000111l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 22,22n|less
./doMafft ROI_regions_tags.txt
./move2LG 
./getAnchors aln/LG26 >LG26_anchors.txt
cp LG25_build.txt LG26_build.txt 
vi LG26_build.txt 
./buildIt LG26_build.txt >LG26.fa
head -2 LG26.fa >LG26_Fa.fa 
tail -2 LG26.fa >LG26_Mo.fa 
cat LG2?_Fa.fa >LG2_Fa.fa
cat LG2?_Mo.fa >LG2_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG2_Fa.fa LG2_Mo.fa >/tmp/m2.paf
minidot /tmp/m2.paf >/tmp/m2.eps
epstopdf /tmp/m2.eps 

cp LG26_build.txt LG27_build.txt 
vi LG27_build.txt 
./buildIt LG27_build.txt >LG27.fa
head -2 LG27.fa >LG27_Fa.fa 
tail -2 LG27.fa >LG27_Mo.fa 
cat LG2?_Fa.fa >LG2_Fa.fa
cat LG2?_Mo.fa >LG2_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG2_Fa.fa LG2_Mo.fa >/tmp/m2.paf
minidot /tmp/m2.paf >/tmp/m2.eps
epstopdf /tmp/m2.eps 

cp LG27_build.txt LG28_build.txt 
vi LG28_build.txt 
./buildIt LG28_build.txt >LG28.fa
head -2 LG28.fa >LG28_Fa.fa 
tail -2 LG28.fa >LG28_Mo.fa 
cat LG2?_Fa.fa >LG2_Fa.fa
cat LG2?_Mo.fa >LG2_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG2_Fa.fa LG2_Mo.fa >/tmp/m2.paf
minidot /tmp/m2.paf >/tmp/m2.eps
epstopdf /tmp/m2.eps 

./getAnchors aln/LG29 >LG29_anchors.txt
perl -ane 'print if $F[19] eq "Mo_000104l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 22,22n|less
./doMafft ROI_regions_tags.txt
./move2LG 
./getAnchors aln/LG29 >LG29_anchors.txt
less LG29_anchors.txt 
cp LG29_anchors.txt LG29_build.txt 
cat LG28_build.txt >>LG29_build.txt 
vi LG29_build.txt 
./buildIt LG29_build.txt >LG29.fa
head -2 LG29.fa >LG29_Fa.fa 
tail -2 LG29.fa >LG29_Mo.fa 
cat LG2?_Fa.fa >LG2_Fa.fa
cat LG2?_Mo.fa >LG2_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG2_Fa.fa LG2_Mo.fa >/tmp/m2.paf
minidot /tmp/m2.paf >/tmp/m2.eps
epstopdf /tmp/m2.eps 

cp LG29_build.txt LG30_build.txt 
vi LG30_build.txt 
./buildIt LG30_build.txt >LG30.fa
head -2 LG30.fa >LG30_Fa.fa 
tail -2 LG30.fa >LG30_Mo.fa 
cat LG3?_Mo.fa >LG3_Mo.fa
cat LG3?_Fa.fa >LG3_Fa.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG3_Fa.fa LG3_Mo.fa >/tmp/m3.paf
minidot /tmp/m3.paf >/tmp/m3.eps
epstopdf /tmp/m3.eps 

cp LG30_build.txt LG31_build.txt 
grep Mo_000024 LG??_build.txt >>LG31_build.txt
vi LG31_build.txt 
./buildIt LG31_build.txt >LG31.fa
head -2 LG31.fa >LG31_Fa.fa 
tail -2 LG31.fa >LG31_Mo.fa 
cat LG3?_Fa.fa >LG3_Fa.fa
cat LG3?_Mo.fa >LG3_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG3_Fa.fa LG3_Mo.fa >/tmp/m3.paf
minidot /tmp/m3.paf >/tmp/m3.eps
epstopdf /tmp/m3.eps 

cp LG31_build.txt LG32_build.txt 
vi LG32_build.txt 
./buildIt LG32_build.txt >LG32.fa
head -2 LG32.fa >LG32_Fa.fa 
tail -2 LG32.fa >LG32_Mo.fa 
cat LG3?_Fa.fa >LG3_Fa.fa
cat LG3?_Mo.fa >LG3_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG3_Fa.fa LG3_Mo.fa >/tmp/m3.paf
minidot /tmp/m3.paf >/tmp/m3.eps
epstopdf /tmp/m3.eps 

./getAnchors aln/LG33 >LG33_anchors.txt
less LG33_anchors.txt 
cp LG33_anchors.txt LG33_build.txt 
cat LG32_build.txt >>LG33_build.txt
vi LG33_build.txt 
./buildIt LG33_build.txt >LG33.fa
head -2 LG33.fa >LG33_Fa.fa 
tail -2 LG33.fa >LG33_Mo.fa 
cat LG3?_Fa.fa >LG3_Fa.fa
cat LG3?_Mo.fa >LG3_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG3_Fa.fa LG3_Mo.fa >/tmp/m3.paf
minidot /tmp/m3.paf >/tmp/m3.eps
epstopdf /tmp/m3.eps 

cp LG33_build.txt LG34_build.txt 
vi LG34_build.txt 
./buildIt LG34_build.txt >LG34.fa
head -2 LG34.fa >LG34_Fa.fa 
tail -2 LG34.fa >LG34_Mo.fa 
cat LG3?_Fa.fa >LG3_Fa.fa
cat LG3?_Mo.fa >LG3_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG3_Fa.fa LG3_Mo.fa >/tmp/m3.paf
minidot /tmp/m3.paf >/tmp/m3.eps
epstopdf /tmp/m3.eps 

cp LG34_build.txt LG35_build.txt 
vi LG35_build.txt 
./buildIt LG35_build.txt >LG35.fa
head -2 LG35.fa >LG35_Fa.fa 
tail -2 LG35.fa >LG35_Mo.fa 
cat LG3?_Fa.fa >LG3_Fa.fa
cat LG3?_Mo.fa >LG3_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG3_Fa.fa LG3_Mo.fa >/tmp/m3.paf
minidot /tmp/m3.paf >/tmp/m3.eps
epstopdf /tmp/m3.eps 

cp LG35_build.txt LG36_build.txt 
vi LG36_build.txt 
./buildIt LG36_build.txt >LG36.fa
head -2 LG36.fa >LG36_Fa.fa 
tail -2 LG36.fa >LG36_Mo.fa 
cat LG3?_Fa.fa >LG3_Fa.fa
cat LG3?_Mo.fa >LG3_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG3_Fa.fa LG3_Mo.fa >/tmp/m3.paf
minidot /tmp/m3.paf >/tmp/m3.eps
epstopdf /tmp/m3.eps 

cp LG36_build.txt LG37_build.txt 
vi LG37_build.txt 
perl -ane 'print if $F[10] eq "Fa_000036l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 13,13n|less
perl -ane 'print if $F[19] eq "Mo_000166l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 22,22n|less
vi LG37_build.txt 
./buildIt LG37_build.txt >LG37.fa
head -2 LG37.fa >LG37_Fa.fa 
tail -2 LG37.fa >LG37_Mo.fa 
cat LG3?_Fa.fa >LG3_Fa.fa
cat LG3?_Mo.fa >LG3_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG3_Fa.fa LG3_Mo.fa >/tmp/m3.paf
minidot /tmp/m3.paf >/tmp/m3.eps
epstopdf /tmp/m3.eps 

./getAnchors aln/LG38 >LG38_anchors.txt
less LG38_anchors.txt 
cp LG37_build.txt LG38_build.txt 
vi LG38_build.txt 
./buildIt LG38_build.txt >LG38.fa
head -2 LG38.fa >LG38_Fa.fa 
tail -2 LG38.fa >LG38_Mo.fa 
cat LG3?_Fa.fa >LG3_Fa.fa
cat LG3?_Mo.fa >LG3_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG3_Fa.fa LG3_Mo.fa >/tmp/m3.paf
minidot /tmp/m3.paf >/tmp/m3.eps
epstopdf /tmp/m3.eps 

cp LG38_build.txt LG39_build.txt 
vi LG39_build.txt 
perl -ane 'print if $F[10] eq "Fa_000006l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 13,13n|less
vi LG39_build.txt 
./buildIt LG39_build.txt >LG39.fa
head -2 LG39.fa >LG39_Fa.fa 
tail -2 LG39.fa >LG39_Mo.fa 
cat LG3?_Fa.fa >LG3_Fa.fa
cat LG3?_Mo.fa >LG3_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG3_Fa.fa LG3_Mo.fa >/tmp/m3.paf
minidot /tmp/m3.paf >/tmp/m3.eps
epstopdf /tmp/m3.eps 

perl -ane 'print if $F[10] eq "Fa_000012l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 13,13n|less
perl -ane 'print if $F[19] eq "Mo_000042l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 22,22n|less
cp LG39_build.txt LG40_build.txt 
vi LG40_build.txt 
./buildIt LG40_build.txt >LG40.fa
head -2 LG40.fa >LG40_Fa.fa 
tail -2 LG40.fa >LG40_Mo.fa 
cat LG4?_Fa.fa >LG4_Fa.fa
cat LG4?_Mo.fa >LG4_Mo.fa
minimap2 -x asm5 --end-bonus 100 -D -t 24 LG4_Fa.fa LG4_Mo.fa >/tmp/m4.paf
minidot /tmp/m4.paf >/tmp/m4.eps
epstopdf /tmp/m4.eps 
```

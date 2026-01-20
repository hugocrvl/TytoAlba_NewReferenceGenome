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

== Mapping of some mRNA of interest and all HiFi reads to the assembled contigs ==

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


hifiasm --version
0.25.0-r726

hifiasm --ont -o hifiasm.Fa/asm -t 40 --primary --n-hap 1 -l 0 --min-hist-cnt 3 ON?.{Father,Unslct}.{other,musat,telo}.fq.gz 2>&1 | tee hifiasm_Fa_log.txt
hifiasm --ont -o hifiasm.Mo/asm -t 40 --primary --n-hap 1 -l 0 --min-hist-cnt 3 ON?.{Mother,Unslct}.{other,musat,telo}.fq.gz 2>&1 | tee hifiasm_Mo_log.txt

sed 's/^>/>ON_Fa_/' hifiasm.Fa/asm.p_ctg.fa >hifiasm_Fa.fa
sed 's/^>/>ON_Mo_/' hifiasm.Mo/asm.p_ctg.fa >hifiasm_Mo.fa

for p1 in Fa Mo; do for p2 in Fa Mo; do echo "$p1 vs $p2"; minimap2 -t 40 -x asm5 -I 12G --secondary=no $p1.0518.p_ctg.fa.gz hifiasm_$p2.fa |perl -ane 'print if ($F[3]-$F[2])>=40000'>$p1.0518.vs.$p2.paf; done; done

./runLastz

./doMafft ROI_regions_tags.txt

./move2LG 

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

 9217  ./getAnchors aln/LG04 >LG04_anchors.txt
 9218  less LG04_anchors.txt 
 9219  cp LG04_anchors.txt LG04_build.txt 
 9220  cat LG03_build.txt >>LG04_build.txt 
 9221  vi LG04_build.txt 
 9222  ./buildIt LG04_build.txt >LG04.fa
 9223  head -2 LG04.fa >LG04_Fa.fa 
 9224  tail -2 LG04.fa >LG04_Mo.fa 
 9225  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG04_Fa.fa LG04_Mo.fa >/tmp/m.paf
 9226  minidot /tmp/m.paf >/tmp/m.eps
 9227  epstopdf /tmp/m.eps 
 9228  cat LG??_Fa.fa >LG_Fa.fa
 9229  cat LG??_Mo.fa >LG_Mo.fa
 9230  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG_Fa.fa LG_Mo.fa >/tmp/m.paf
 9231  minidot /tmp/m.paf >/tmp/m.eps
 9232  epstopdf /tmp/m.eps 

 9233  ./getAnchors aln/LG05 >LG05_anchors.txt
 9234  less LG05_anchors.txt 
 9235  perl -ane 'print if $F[10] eq "Mo_000019l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 13,13n|less
 9236  head -1 /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt
 9237  perl -ane 'print if $F[19] eq "Mo_000019l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 22,22n|less
 9238  vi doMafft 
 9239  ./doMafft ROI_regions_tags.txt
 9240  rm aln/Fa_000005l_12252186_12287599_all.*
 9241  vi ROI_regions_tags.txt 
 9242  ./doMafft ROI_regions_tags.txt
 9243  rm aln/Mo_000019l_31567601_31606358_*
 9245  ./move2LG 
 9246  ./getAnchors aln/LG05 >LG05_anchors.txt
 9247  less LG05_anchors.txt 
 9248  cp LG05_anchors.txt LG05_build.txt 
 9249  cat LG05_build.txt >>LG05_build.txt 
 9250  cat LG04_build.txt >>LG05_build.txt 
 9251  vi LG05_build.txt 
 9252  ./buildIt LG05_build.txt >LG05.fa
 9253  head -2 LG05.fa >LG05_Fa.fa 
 9254  tail -2 LG05.fa >LG05_Mo.fa 
 9255  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG05_Fa.fa LG05_Mo.fa >/tmp/m.paf
 9256  minidot /tmp/m.paf >/tmp/m.eps
 9257  epstopdf /tmp/m.eps 

 9258  ./getAnchors aln/LG06 >LG06_anchors.txt
 9259  less LG06_anchors.txt 
 9260  perl -ane 'print if $F[19] eq "Mo_000005l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 22,22n|less
 9261  ./doMafft ROI_regions_tags.txt
 9262  ./move2LG 
 9263  ./getAnchors aln/LG06 >LG06_anchors.txt
 9264  less LG06_anchors.txt 
 9265  cp LG06_anchors.txt LG06_build.txt 
 9266  cat LG05_build.txt >>LG06_build.txt 
 9267  vi LG06_build.txt 
 9268  grep Fa_000132l /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt|head|less
 9269  vi LG06_build.txt 
 9270  ./buildIt LG06_build.txt >LG06.fa
 9271  head -2 LG06.fa >LG06_Fa.fa 
 9272  tail -2 LG06.fa >LG06_Mo.fa 
 9273  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG06_Fa.fa LG06_Mo.fa >/tmp/m.paf
 9274  minidot /tmp/m.paf >/tmp/m.eps
 9275  epstopdf /tmp/m.eps 
 9276  cat LG??_Fa.fa >LG_Fa.fa
 9277  cat LG??_Mo.fa >LG_Mo.fa
 9278  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG_Fa.fa LG_Mo.fa >/tmp/m.paf
 9279  minidot /tmp/m.paf >/tmp/m.eps
 9280  epstopdf /tmp/m.eps 

 9281  cp LG06_build.txt LG07_build.txt 
 9282  vi LG07_build.txt 
 9283  grep Mo_000049l /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt|head -2 >>LG07_build.txt 
 9284  grep Fa_000015l /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt|head -2 >>LG07_build.txt 
 9285  vi LG07_build.txt 
 9286  ./buildIt LG07_build.txt >LG07.fa
 9287  head -2 LG07.fa >LG07_Fa.fa 
 9288  tail -2 LG07.fa >LG07_Mo.fa 
 9289  cat LG??_Fa.fa >LG_Fa.fa
 9290  cat LG??_Mo.fa >LG_Mo.fa
 9291  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG_Fa.fa LG_Mo.fa >/tmp/m.paf
 9292  minidot /tmp/m.paf >/tmp/m.eps
 9293  epstopdf /tmp/m.eps 

 9294  cp LG07_build.txt LG08_build.txt 
 9295  grep Mo_000017l /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt|head -2 >>LG08_build.txt 
 9296  vi LG08_build.txt 
 9297  ./buildIt LG08_build.txt >LG08.fa
 9298  head -2 LG08.fa >LG08_Fa.fa 
 9299  tail -2 LG08.fa >LG08_Mo.fa 
 9300  cat LG??_Fa.fa >LG_Fa.fa
 9301  cat LG??_Mo.fa >LG_Mo.fa
 9302  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG_Fa.fa LG_Mo.fa >/tmp/m.paf
 9303  minidot /tmp/m.paf >/tmp/m.eps
 9304  epstopdf /tmp/m.eps 

 9305  ./getAnchors aln/LG09 >LG09_anchors.txt
 9306  less LG09_anchors.txt 
 9307  perl -ane 'print if $F[19] eq "Mo_000075l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 22,22n|less
 9308  ./doMafft ROI_regions_tags.txt
 9309  ./move2LG 
 9310  ./getAnchors aln/LG09 >LG09_anchors.txt
 9311  less LG09_anchors.txt 
 9312  cp LG09_anchors.txt LG09_build.txt 
 9313  cat LG08_build.txt >>LG09_build.txt 
 9314  vi LG09_build.txt 
 9315  ./buildIt LG09_build.txt >LG09.fa
 9316  head -2 LG09.fa >LG09_Fa.fa 
 9317  tail -2 LG09.fa >LG09_Mo.fa 
 9318  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG09_Fa.fa LG09_Mo.fa >/tmp/m.paf
 9319  minidot /tmp/m.paf >/tmp/m.eps
 9320  epstopdf /tmp/m.eps 

 9321  ./getAnchors aln/LG10 >LG10_anchors.txt
 9322  less LG10_anchors.txt 
 9323  perl -ane 'print if $F[10] eq "Fa_000105l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 13,13n|less
 9324  ./doMafft ROI_regions_tags.txt
 9325  ./move2LG 
 9326  ./getAnchors aln/LG10 >LG10_anchors.txt
 9327  less LG10_anchors.txt 
 9328  ls aln/LG10/
 9329  less aln/LG10/Fa_000105l_22927074_22938956_log.txt 
 9330  cp LG09_build.txt LG10_build.txt 
 9331  vi LG10_build.txt 
 9332  cat LG10_anchors.txt >>LG10_build.txt 
 9333  vi LG10_build.txt 
 9334  ./buildIt LG10_build.txt >LG10.fa
 9335  head -2 LG10.fa >LG10_Fa.fa 
 9336  tail -2 LG10.fa >LG10_Mo.fa 
 9337  cat LG??_Fa.fa >LG_Fa.fa
 9338  cat LG??_Mo.fa >LG_Mo.fa
 9339  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG_Fa.fa LG_Mo.fa >/tmp/m.paf
 9340  minidot /tmp/m.paf >/tmp/m.eps
 9341  epstopdf /tmp/m.eps 
 9342  vi LG10_build.txt 
 9343  ./buildIt LG10_build.txt >LG10.fa
 9344  head -2 LG10.fa >LG10_Fa.fa 
 9345  tail -2 LG10.fa >LG10_Mo.fa 
 9346  cat LG??_Fa.fa >LG_Fa.fa
 9347  cat LG??_Mo.fa >LG_Mo.fa
 9348  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG_Fa.fa LG_Mo.fa >/tmp/m.paf
 9349  minidot /tmp/m.paf >/tmp/m.eps
 9350  epstopdf /tmp/m.eps 

 9362  ./getAnchors aln/LG11 >LG11_anchors.txt
 9363  less LG11_anchors.txt 
 9364  cp LG11_anchors.txt LG11_build.txt 
 9365  cat LG10_build.txt >>LG11_build.txt 
 9366  vi LG11_build.txt 
 9367  ./buildIt LG11_build.txt >LG11.fa
 9368  head -2 LG11.fa >LG11_Fa.fa 
 9369  tail -2 LG11.fa >LG11_Mo.fa 
 9370  cat LG??_Fa.fa >LG_Fa.fa
 9371  cat LG??_Mo.fa >LG_Mo.fa
 9372  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG_Fa.fa LG_Mo.fa >/tmp/m.paf
 9373  minidot /tmp/m.paf >/tmp/m.eps
 9374  epstopdf /tmp/m.eps 

 9381  ./getAnchors aln/LG12 >LG12_anchors.txt
 9382  less LG12_anchors.txt 
 9383  cp LG12_anchors.txt LG12_build.txt 
 9384  cat LG11_build.txt >>LG12_build.txt 
 9385  vi LG12_build.txt 
 9386  ./buildIt LG12_build.txt >LG12.fa
 9387  cat LG??_Fa.fa >LG_Fa.fa
 9388  cat LG??_Mo.fa >LG_Mo.fa
 9389  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG_Fa.fa LG_Mo.fa >/tmp/m.paf
 9390  minidot /tmp/m.paf >/tmp/m.eps
 9391  epstopdf /tmp/m.eps 
 9392  head -2 LG12.fa >LG12_Fa.fa 
 9393  tail -2 LG12.fa >LG12_Mo.fa 
 9394  cat LG??_Mo.fa >LG_Mo.fa
 9395  cat LG??_Fa.fa >LG_Fa.fa
 9396  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG_Fa.fa LG_Mo.fa >/tmp/m.paf
 9397  minidot /tmp/m.paf >/tmp/m.eps
 9398  epstopdf /tmp/m.eps 

 9399  cp LG07_build.txt LG13_build.txt 
 9400  vi LG13_build.txt 
 9401  ./buildIt LG13_build.txt >LG13.fa
 9402  head -2 LG13.fa >LG13_Fa.fa 
 9403  tail -2 LG13.fa >LG13_Mo.fa 
 9404  cat LG??_Fa.fa >LG_Fa.fa
 9405  cat LG??_Mo.fa >LG_Mo.fa
 9406  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG_Fa.fa LG_Mo.fa >/tmp/m.paf
 9407  minidot /tmp/m.paf >/tmp/m.eps
 9408  epstopdf /tmp/m.eps 
 9409  cat LG0?_Fa.fa >LG0_Fa.fa
 9410  cat LG0?_Mo.fa >LG0_Mo.fa
 9411  cat LG1?_Fa.fa >LG1_Fa.fa
 9412  cat LG1?_Mo.fa >LG1_Mo.fa
 9413  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG0_Fa.fa LG0_Mo.fa >/tmp/m0.paf
 9414  minidot /tmp/m0.paf >/tmp/m0.eps
 9415  epstopdf /tmp/m0.eps 
 9416  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG1_Fa.fa LG1_Mo.fa >/tmp/m1.paf
 9417  minidot /tmp/m1.paf >/tmp/m1.eps
 9418  epstopdf /tmp/m1.eps 

 9419  cp LG13_build.txt LG14_build.txt 
 9420  vi LG14_build.txt 
 9421  ./buildIt LG14_build.txt >LG14.fa
 9422  head -2 LG14.fa >LG14_Fa.fa 
 9423  tail -2 LG14.fa >LG14_Mo.fa 
 9424  cat LG1?_Fa.fa >LG1_Fa.fa
 9425  cat LG1?_Mo.fa >LG1_Mo.fa
 9426  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG1_Fa.fa LG1_Mo.fa >/tmp/m1.paf
 9427  minidot /tmp/m1.paf >/tmp/m1.eps
 9428  epstopdf /tmp/m1.eps 

 9434  ./getAnchors aln/LG15 >LG15_anchors.txt
 9435  less LG15_anchors.txt 
 9436  cp LG15_anchors.txt LG15_build.txt 
 9437  cat LG14_build.txt >>LG15_build.txt 
 9438  vi LG15_build.txt 
 9439  ./buildIt LG15_build.txt >LG15.fa
 9440  cat LG1?_Fa.fa >LG1_Fa.fa
 9441  cat LG1?_Mo.fa >LG1_Mo.fa
 9442  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG1_Fa.fa LG1_Mo.fa >/tmp/m1.paf
 9443  minidot /tmp/m1.paf >/tmp/m1.eps
 9444  epstopdf /tmp/m1.eps 
 9445  head -2 LG15.fa >LG15_Fa.fa 
 9446  tail -2 LG15.fa >LG15_Mo.fa 
 9447  cat LG1?_Fa.fa >LG1_Fa.fa
 9448  cat LG1?_Mo.fa >LG1_Mo.fa
 9449  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG1_Fa.fa LG1_Mo.fa >/tmp/m1.paf
 9450  minidot /tmp/m1.paf >/tmp/m1.eps
 9451  epstopdf /tmp/m1.eps 

 9452  ./getAnchors aln/LG16 >LG16_anchors.txt
 9453  perl -ane 'print if $F[19] eq "Mo_000001l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 22,22n|less
 9454  ./doMafft ROI_regions_tags.txt
 9455  ./move2LG 
 9456  ./getAnchors aln/LG16 >LG16_anchors.txt
 9457  less LG16_anchors.txt 
 9458  cp LG16_anchors.txt LG16_build.txt 
 9459  cat LG15_build.txt >>LG16_build.txt 
 9460  vi LG16_build.txt 
 9461  ./buildIt LG16_build.txt >LG16.fa
 9462  head -2 LG65.fa >LG16_Fa.fa 
 9463  head -2 LG16.fa >LG16_Fa.fa 
 9464  tail -2 LG16.fa >LG16_Mo.fa 
 9465  cat LG1?_Fa.fa >LG1_Fa.fa
 9466  cat LG1?_Mo.fa >LG1_Mo.fa
 9467  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG1_Fa.fa LG1_Mo.fa >/tmp/m1.paf
 9468  minidot /tmp/m1.paf >/tmp/m1.eps
 9469  epstopdf /tmp/m1.eps 
 9470  vi /scratch/tyto/plot_cM_pos_11_20.R 
 9471  vi LG16_build.txt 
 9472  ./buildIt LG16_build.txt >LG16.fa
 9473  head -2 LG16.fa >LG16_Fa.fa 
 9474  tail -2 LG16.fa >LG16_Mo.fa 
 9475  cat LG1?_Fa.fa >LG1_Fa.fa
 9476  cat LG1?_Mo.fa >LG1_Mo.fa
 9477  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG1_Fa.fa LG1_Mo.fa >/tmp/m1.paf
 9478  minidot /tmp/m1.paf >/tmp/m1.eps
 9479  epstopdf /tmp/m1.eps 

 9480  ./getAnchors aln/LG17 >LG17_anchors.txt
 9481  less LG17_anchors.txt 
 9482  perl -ane 'print if $F[10] eq "Fa_000057l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 13,13n|less
 9483  ./doMafft ROI_regions_tags.txt
 9484  ./move2LG 
 9485  ./getAnchors aln/LG17 >LG17_anchors.txt
 9486  cp LG17_anchors.txt LG17_build.txt 
 9487  cat LG16_build.txt >>LG17_build.txt 
 9488  vi LG17_build.txt 
 9489  ./buildIt LG17_build.txt >LG17.fa
 9490  head -2 LG17.fa >LG17_Fa.fa 
 9491  tail -2 LG17.fa >LG17_Mo.fa 
 9492  cat LG1?_Fa.fa >LG1_Fa.fa
 9493  cat LG1?_Mo.fa >LG1_Mo.fa
 9494  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG1_Fa.fa LG1_Mo.fa >/tmp/m1.paf
 9495  minidot /tmp/m1.paf >/tmp/m1.eps
 9496  epstopdf /tmp/m1.eps 
 9497  vi /scratch/tyto/plot_cM_pos_11_20.R 
 9498  vi LG17_build.txt 
 9499  ./buildIt LG17_build.txt >LG17.fa
 9500  head -2 LG17.fa >LG17_Fa.fa 
 9501  tail -2 LG17.fa >LG17_Mo.fa 
 9502  cat LG1?_Fa.fa >LG1_Fa.fa
 9503  cat LG1?_Mo.fa >LG1_Mo.fa
 9504  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG1_Fa.fa LG1_Mo.fa >/tmp/m1.paf
 9505  minidot /tmp/m1.paf >/tmp/m1.eps
 9506  epstopdf /tmp/m1.eps 

 9507  cp LG17_build.txt LG18_build.txt 
 9508  vi LG18_build.txt 
 9509  perl -ane 'print if $F[19] eq "Mo_000003l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 22,22n|less
 9510  ls -R|grep Fa_000091
 9511  perl -ane 'print if $F[19] eq "Mo_000003l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 22,22n|less
 9512  vi LG18_build.txt 
 9513  ./buildIt LG18_build.txt >LG18.fa
 9514  head -2 LG18.fa >LG18_Fa.fa 
 9515  tail -2 LG18.fa >LG18_Mo.fa 
 9516  cat LG1?_Fa.fa >LG1_Fa.fa
 9517  cat LG1?_Mo.fa >LG1_Mo.fa
 9518  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG1_Fa.fa LG1_Mo.fa >/tmp/m1.paf
 9519  minidot /tmp/m1.paf >/tmp/m1.eps
 9520  epstopdf /tmp/m1.eps 

 9521  cp LG18_build.txt LG19_build.txt 
 9522  vi LG19_build.txt 
 9523  ./buildIt LG19_build.txt >LG19.fa
 9524  head -2 LG19.fa >LG19_Fa.fa 
 9525  tail -2 LG19.fa >LG19_Mo.fa 
 9526  cat LG1?_Fa.fa >LG1_Fa.fa
 9527  cat LG1?_Mo.fa >LG1_Mo.fa
 9528  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG1_Fa.fa LG1_Mo.fa >/tmp/m1.paf
 9529  minidot /tmp/m1.paf >/tmp/m1.eps
 9530  epstopdf /tmp/m1.eps 

 9531  cp LG19_build.txt LG20_build.txt 
 9532  vi LG20_build.txt 
 9533  ./buildIt LG20_build.txt >LG20.fa
 9534  head -2 LG20.fa >LG20_Fa.fa 
 9535  tail -2 LG20.fa >LG20_Mo.fa 
 9536  cat LG2?_Fa.fa >LG2_Fa.fa
 9537  cat LG2?_Mo.fa >LG2_Mo.fa
 9538  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG2_Fa.fa LG2_Mo.fa >/tmp/m2.paf
 9539  minidot /tmp/m1.paf >/tmp/m2.eps
 9540  minidot /tmp/m2.paf >/tmp/m2.eps
 9541  epstopdf /tmp/m2.eps 

 9542  cp LG20_build.txt LG21_build.txt 
 9543  vi LG21_build.txt 
 9544  perl -ane 'print if $F[10] eq "Fa_000053l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 13,13n|less
 9545  vi LG21_build.txt 
 9546  ./buildIt LG21_build.txt >LG21.fa
 9547  head -2 LG21.fa >LG21_Fa.fa 
 9548  tail -2 LG21.fa >LG21_Mo.fa 
 9549  cat LG2?_Fa.fa >LG2_Fa.fa
 9550  cat LG2?_Mo.fa >LG2_Mo.fa
 9551  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG2_Fa.fa LG2_Mo.fa >/tmp/m2.paf
 9552  minidot /tmp/m2.paf >/tmp/m2.eps
 9553  epstopdf /tmp/m2.eps 

 9554  cp LG21_build.txt LG22_build.txt 
 9555  vi LG22_build.txt 
 9556  ./buildIt LG22_build.txt >LG22.fa
 9557  head -2 LG22.fa >LG22_Fa.fa 
 9558  tail -2 LG22.fa >LG22_Mo.fa 
 9559  cat LG2?_Fa.fa >LG2_Fa.fa
 9560  cat LG2?_Mo.fa >LG2_Mo.fa
 9561  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG2_Fa.fa LG2_Mo.fa >/tmp/m2.paf
 9562  minidot /tmp/m2.paf >/tmp/m2.eps
 9563  epstopdf /tmp/m2.eps 

 9564  ./getAnchors aln/LG23 >LG23_anchors.txt
 9565  perl -ane 'print if $F[10] eq "Fa_000167l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 13,13n|less
 9566  ./doMafft ROI_regions_tags.txt
 9567  ./getAnchors aln/LG23 >LG23_anchors.txt
 9568  ./move2LG 
 9569  ./getAnchors aln/LG23 >LG23_anchors.txt
 9570  less LG23_anchors.txt 
 9571  cp LG22_build.txt LG23_build.txt 
 9572  vi LG23_build.txt 
 9573  ./buildIt LG23_build.txt >LG23.fa
 9574  head -2 LG23.fa >LG23_Fa.fa 
 9575  tail -2 LG23.fa >LG23_Mo.fa 
 9576  cat LG2?_Fa.fa >LG2_Fa.fa
 9577  cat LG2?_Mo.fa >LG2_Mo.fa
 9578  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG2_Fa.fa LG2_Mo.fa >/tmp/m2.paf
 9579  minidot /tmp/m2.paf >/tmp/m2.eps
 9580  epstopdf /tmp/m2.eps 

 9581  ./getAnchors aln/LG24 >LG24_anchors.txt
 9582  less LG24_anchors.txt 
 9583  perl -ane 'print if $F[19] eq "Mo_000016l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 22,22n|less
 9584  ./doMafft ROI_regions_tags.txt
 9585  ./move2LG 
 9586  ./getAnchors aln/LG24 >LG24_anchors.txt
 9587  less LG24_anchors.txt 
 9588  cp LG24_anchors.txt LG24_build.txt 
 9589  cat LG23_build.txt >>LG24_build.txt 
 9590  vi LG24_build.txt 
 9591  ./buildIt LG24_build.txt >LG24.fa
 9592  head -2 LG24.fa >LG24_Fa.fa 
 9593  tail -2 LG24.fa >LG24_Mo.fa 
 9594  cat LG2?_Fa.fa >LG2_Fa.fa
 9595  cat LG2?_Mo.fa >LG2_Mo.fa
 9596  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG2_Fa.fa LG2_Mo.fa >/tmp/m2.paf
 9597  minidot /tmp/m2.paf >/tmp/m2.eps
 9598  epstopdf /tmp/m2.eps 

 9599  cp LG24_build.txt LG25_build.txt 
 9600  vi LG25_build.txt 
 9601  ./buildIt LG25_build.txt >LG25.fa
 9602  head -2 LG25.fa >LG25_Fa.fa 
 9603  tail -2 LG25.fa >LG25_Mo.fa 
 9604  cat LG2?_Fa.fa >LG2_Fa.fa
 9605  cat LG2?_Mo.fa >LG2_Mo.fa
 9606  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG2_Fa.fa LG2_Mo.fa >/tmp/m2.paf
 9607  minidot /tmp/m2.paf >/tmp/m2.eps
 9608  epstopdf /tmp/m2.eps 

 9609  ./getAnchors aln/LG26 >LG26_anchors.txt
 9610  perl -ane 'print if $F[10] eq "Fa_000065l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 13,13n|less
 9611  perl -ane 'print if $F[19] eq "Mo_000111l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 22,22n|less
 9612  ./doMafft ROI_regions_tags.txt
 9613  ./move2LG 
 9614  ./getAnchors aln/LG26 >LG26_anchors.txt
 9615  cp LG25_build.txt LG26_build.txt 
 9616  vi LG26_build.txt 
 9617  ./buildIt LG26_build.txt >LG26.fa
 9618  head -2 LG26.fa >LG26_Fa.fa 
 9619  tail -2 LG26.fa >LG26_Mo.fa 
 9620  cat LG2?_Fa.fa >LG2_Fa.fa
 9621  cat LG2?_Mo.fa >LG2_Mo.fa
 9622  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG2_Fa.fa LG2_Mo.fa >/tmp/m2.paf
 9623  minidot /tmp/m2.paf >/tmp/m2.eps
 9624  epstopdf /tmp/m2.eps 

 9625  cp LG26_build.txt LG27_build.txt 
 9626  vi LG27_build.txt 
 9627  ./buildIt LG27_build.txt >LG27.fa
 9628  head -2 LG27.fa >LG27_Fa.fa 
 9629  tail -2 LG27.fa >LG27_Mo.fa 
 9630  cat LG2?_Fa.fa >LG2_Fa.fa
 9631  cat LG2?_Mo.fa >LG2_Mo.fa
 9632  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG2_Fa.fa LG2_Mo.fa >/tmp/m2.paf
 9633  minidot /tmp/m2.paf >/tmp/m2.eps
 9634  epstopdf /tmp/m2.eps 

 9635  cp LG27_build.txt LG28_build.txt 
 9636  vi LG28_build.txt 
 9637  ./buildIt LG28_build.txt >LG28.fa
 9638  head -2 LG28.fa >LG28_Fa.fa 
 9639  tail -2 LG28.fa >LG28_Mo.fa 
 9640  cat LG2?_Fa.fa >LG2_Fa.fa
 9641  cat LG2?_Mo.fa >LG2_Mo.fa
 9642  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG2_Fa.fa LG2_Mo.fa >/tmp/m2.paf
 9643  minidot /tmp/m2.paf >/tmp/m2.eps
 9644  epstopdf /tmp/m2.eps 

 9645  ./getAnchors aln/LG29 >LG29_anchors.txt
 9646  perl -ane 'print if $F[19] eq "Mo_000104l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 22,22n|less
 9647  ./doMafft ROI_regions_tags.txt
 9648  ./move2LG 
 9649  ./getAnchors aln/LG29 >LG29_anchors.txt
 9650  less LG29_anchors.txt 
 9651  cp LG29_anchors.txt LG29_build.txt 
 9652  cat LG28_build.txt >>LG29_build.txt 
 9653  vi LG29_build.txt 
 9654  ./buildIt LG29_build.txt >LG29.fa
 9655  head -2 LG29.fa >LG29_Fa.fa 
 9656  tail -2 LG29.fa >LG29_Mo.fa 
 9657  cat LG2?_Fa.fa >LG2_Fa.fa
 9658  cat LG2?_Mo.fa >LG2_Mo.fa
 9659  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG2_Fa.fa LG2_Mo.fa >/tmp/m2.paf
 9660  minidot /tmp/m2.paf >/tmp/m2.eps
 9661  epstopdf /tmp/m2.eps 

 9663  cp LG29_build.txt LG30_build.txt 
 9664  vi LG30_build.txt 
 9665  ./buildIt LG30_build.txt >LG30.fa
 9666  head -2 LG30.fa >LG30_Fa.fa 
 9667  tail -2 LG30.fa >LG30_Mo.fa 
 9668  cat LG3?_Mo.fa >LG3_Mo.fa
 9669  cat LG3?_Fa.fa >LG3_Fa.fa
 9670  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG3_Fa.fa LG3_Mo.fa >/tmp/m3.paf
 9671  minidot /tmp/m3.paf >/tmp/m3.eps
 9672  epstopdf /tmp/m3.eps 

 9674  cp LG30_build.txt LG31_build.txt 
 9675  grep Mo_000024 LG??_build.txt >>LG31_build.txt
 9676  vi LG31_build.txt 
 9677  ./buildIt LG31_build.txt >LG31.fa
 9678  head -2 LG31.fa >LG31_Fa.fa 
 9679  tail -2 LG31.fa >LG31_Mo.fa 
 9680  cat LG3?_Fa.fa >LG3_Fa.fa
 9681  cat LG3?_Mo.fa >LG3_Mo.fa
 9682  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG3_Fa.fa LG3_Mo.fa >/tmp/m3.paf
 9683  minidot /tmp/m3.paf >/tmp/m3.eps
 9684  epstopdf /tmp/m3.eps 

 9685  cp LG31_build.txt LG32_build.txt 
 9686  vi LG32_build.txt 
 9687  ./buildIt LG32_build.txt >LG32.fa
 9688  head -2 LG32.fa >LG32_Fa.fa 
 9689  tail -2 LG32.fa >LG32_Mo.fa 
 9690  cat LG3?_Fa.fa >LG3_Fa.fa
 9691  cat LG3?_Mo.fa >LG3_Mo.fa
 9692  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG3_Fa.fa LG3_Mo.fa >/tmp/m3.paf
 9693  minidot /tmp/m3.paf >/tmp/m3.eps
 9694  epstopdf /tmp/m3.eps 

 9695  ./getAnchors aln/LG33 >LG33_anchors.txt
 9696  less LG33_anchors.txt 
 9697  cp LG33_anchors.txt LG33_build.txt 
 9698  cat LG32_build.txt >>LG33_build.txt
 9699  vi LG33_build.txt 
 9700  ./buildIt LG33_build.txt >LG33.fa
 9701  head -2 LG33.fa >LG33_Fa.fa 
 9702  tail -2 LG33.fa >LG33_Mo.fa 
 9703  cat LG3?_Fa.fa >LG3_Fa.fa
 9704  cat LG3?_Mo.fa >LG3_Mo.fa
 9705  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG3_Fa.fa LG3_Mo.fa >/tmp/m3.paf
 9706  minidot /tmp/m3.paf >/tmp/m3.eps
 9707  epstopdf /tmp/m3.eps 

 9708  cp LG33_build.txt LG34_build.txt 
 9709  vi LG34_build.txt 
 9710  ./buildIt LG34_build.txt >LG34.fa
 9711  head -2 LG34.fa >LG34_Fa.fa 
 9712  tail -2 LG34.fa >LG34_Mo.fa 
 9713  cat LG3?_Fa.fa >LG3_Fa.fa
 9714  cat LG3?_Mo.fa >LG3_Mo.fa
 9715  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG3_Fa.fa LG3_Mo.fa >/tmp/m3.paf
 9716  minidot /tmp/m3.paf >/tmp/m3.eps
 9717  epstopdf /tmp/m3.eps 

 9718  cp LG34_build.txt LG35_build.txt 
 9719  vi LG35_build.txt 
 9720  ./buildIt LG35_build.txt >LG35.fa
 9721  head -2 LG35.fa >LG35_Fa.fa 
 9722  tail -2 LG35.fa >LG35_Mo.fa 
 9723  cat LG3?_Fa.fa >LG3_Fa.fa
 9724  cat LG3?_Mo.fa >LG3_Mo.fa
 9725  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG3_Fa.fa LG3_Mo.fa >/tmp/m3.paf
 9726  minidot /tmp/m3.paf >/tmp/m3.eps
 9727  epstopdf /tmp/m3.eps 

 9728  cp LG35_build.txt LG36_build.txt 
 9729  vi LG36_build.txt 
 9730  ./buildIt LG36_build.txt >LG36.fa
 9731  head -2 LG36.fa >LG36_Fa.fa 
 9732  tail -2 LG36.fa >LG36_Mo.fa 
 9733  cat LG3?_Fa.fa >LG3_Fa.fa
 9734  cat LG3?_Mo.fa >LG3_Mo.fa
 9735  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG3_Fa.fa LG3_Mo.fa >/tmp/m3.paf
 9736  minidot /tmp/m3.paf >/tmp/m3.eps
 9737  epstopdf /tmp/m3.eps 

 9738  cp LG36_build.txt LG37_build.txt 
 9739  vi LG37_build.txt 
 9740  perl -ane 'print if $F[10] eq "Fa_000036l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 13,13n|less
 9741  perl -ane 'print if $F[19] eq "Mo_000166l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 22,22n|less
 9742  vi LG37_build.txt 
 9743  ./buildIt LG37_build.txt >LG37.fa
 9744  head -2 LG37.fa >LG37_Fa.fa 
 9745  tail -2 LG37.fa >LG37_Mo.fa 
 9746  cat LG3?_Fa.fa >LG3_Fa.fa
 9747  cat LG3?_Mo.fa >LG3_Mo.fa
 9748  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG3_Fa.fa LG3_Mo.fa >/tmp/m3.paf
 9749  minidot /tmp/m3.paf >/tmp/m3.eps
 9750  epstopdf /tmp/m3.eps 

 9751  ./getAnchors aln/LG38 >LG38_anchors.txt
 9752  less LG38_anchors.txt 
 9753  cp LG37_build.txt LG38_build.txt 
 9754  vi LG38_build.txt 
 9755  ./buildIt LG38_build.txt >LG38.fa
 9756  head -2 LG38.fa >LG38_Fa.fa 
 9757  tail -2 LG38.fa >LG38_Mo.fa 
 9758  cat LG3?_Fa.fa >LG3_Fa.fa
 9759  cat LG3?_Mo.fa >LG3_Mo.fa
 9760  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG3_Fa.fa LG3_Mo.fa >/tmp/m3.paf
 9761  minidot /tmp/m3.paf >/tmp/m3.eps
 9762  epstopdf /tmp/m3.eps 

 9763  cp LG38_build.txt LG39_build.txt 
 9764  vi LG39_build.txt 
 9765  perl -ane 'print if $F[10] eq "Fa_000006l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 13,13n|less
 9766  vi LG39_build.txt 
 9767  ./buildIt LG39_build.txt >LG39.fa
 9768  head -2 LG39.fa >LG39_Fa.fa 
 9769  tail -2 LG39.fa >LG39_Mo.fa 
 9770  cat LG3?_Fa.fa >LG3_Fa.fa
 9771  cat LG3?_Mo.fa >LG3_Mo.fa
 9772  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG3_Fa.fa LG3_Mo.fa >/tmp/m3.paf
 9773  minidot /tmp/m3.paf >/tmp/m3.eps
 9774  epstopdf /tmp/m3.eps 

 9775  perl -ane 'print if $F[10] eq "Fa_000012l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 13,13n|less
 9776  perl -ane 'print if $F[19] eq "Mo_000042l"' /scratch/tyto/m64_all_tab_map_MoFa_contig_0518.txt |sort -k 22,22n|less
 9777  cp LG39_build.txt LG40_build.txt 
 9778  vi LG40_build.txt 
 9779  ./buildIt LG40_build.txt >LG40.fa
 9780  head -2 LG40.fa >LG40_Fa.fa 
 9781  tail -2 LG40.fa >LG40_Mo.fa 
 9782  cat LG4?_Fa.fa >LG4_Fa.fa
 9783  cat LG4?_Mo.fa >LG4_Mo.fa
 9784  minimap2 -x asm5 --end-bonus 100 -D -t 24 LG4_Fa.fa LG4_Mo.fa >/tmp/m4.paf
 9785  minidot /tmp/m4.paf >/tmp/m4.eps
 9786  epstopdf /tmp/m4.eps 

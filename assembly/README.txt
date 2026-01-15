Singularity> hifiasm  --version
0.19.8-r603

hifiasm -o /data/chris/tyto/20240518Fa/asm -t 72 --primary -D 10.0 -N 1000 --max-kocc 20000 --n-hap 1 --ul SOL000x_SUP_20k_cutdup.fq.gz,LR_40k.fq.gz /data/chris/tyto/20240518Fa/{other,musat,telomers}.fq.gz /data/chris/tyto/20240518Fa/Unslct.{other,musat,telomers}.fq.gz 2>&1 | tee 20240518Fa_log.txt
hifiasm -o /data/chris/tyto/20240518Mo/asm -t 72 --primary -D 10.0 -N 1000 --max-kocc 20000 --n-hap 1 --ul SOL000x_SUP_20k_cutdup.fq.gz,LR_40k.fq.gz /data/chris/tyto/20240518Mo/{other,musat,telomers}.fq.gz /data/chris/tyto/20240518Fa/Unslct.{other,musat,telomers}.fq.gz 2>&1 | tee 20240518Mo_log.txt

pushd /data/chris/tyto/20240518Fa/
minimap2 -x splice:hq asm.p_ctg.fa ~/tyto/SeqGeneForChristian.fasta >Fa_mRNA.paf
cd ../20240518Mo/
minimap2 -x splice:hq asm.p_ctg.fa ~/tyto/SeqGeneForChristian.fasta >Mo_mRNA.paf

./doMapMoFa 

sort -m -k 6,6 -k 8,8n -S 10% /data/chris/tyto/20240518Fa/m64???_*.paf >m64_all_Fa.paf &
sort -m -k 6,6 -k 8,8n -S 10% /data/chris/tyto/20240518Mo/m64???_*.paf >m64_all_Mo.paf &

./genCoverage m64_all_Fa.paf >m64_all_Fa_cov.txt &
./genCoverage m64_all_Mo.paf >m64_all_Mo_cov.txt &

./analyzeCov >zero_cov_regions.txt

./grabCandidates zero_cov_regions.txt >zero_cov_regions_tags.txt


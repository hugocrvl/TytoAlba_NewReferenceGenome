# Miniprot & Sankey analysis

Well-annotated bird species proteins (zebra finch and chicken) were aligned against the new barn owl haplotype assemblies using Miniprot v0.18-r281. Following the approach used in MicroFinder, we filtered alignments to retain only high-quality matches with ≥70% sequence identity and an additional threshold of alignments spanning ≥90% of the target sequence length. Chromosomal location of the protein was parsed from the feature tables, obtained from the same download directories as the amino acid FASTA files. The syntenic matches across species' chromosomes were then visualized using Sankey diagrams from the networkD3 library in R. The next sectinos describe these steps with further detail and link to the corresponding scripts.

## Miniprot analysis

### Input files

Protein FASTA links, with last accession date:

- Chicken (19.12.25): <https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_protein.faa.gz>
- Zebra finch (19.12.25): <https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/771/995/GCF_048771995.1_bTaeGut7.mat/GCF_048771995.1_bTaeGut7.mat_protein.faa.gz>

Protein entries per fasta file:

|Species    | n   |
|-----------|-----|
|Chicken    |18372|
|Zebra finch|16489|

After running miniprot, we combine the results to feature annotation also available on NCBI's page. Download links & last accession date:

- Chicken (19.12.25): <https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_feature_table.txt.gz>
- Zebra finch (19.12.25): <https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/771/995/GCF_048771995.1_bTaeGut7.mat/GCF_048771995.1_bTaeGut7.mat_feature_table.txt.gz>

### Running Miniprot and processing results

```bash
chicken_proteins="GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_protein.faa.gz"
zfinch_proteins="GCF_048771995.1_bTaeGut7.mat_protein.faa.gz"
barn_owl_fa="all_Fa.fa"
barn_owl_mo="all_Mo.fa"

miniprot -t16 $barn_owl_mo $chicken_proteins > barnowlMO_chicken.paf
miniprot -t16 $barn_owl_fa $chicken_proteins > barnowl_chicken.paf

miniprot -t16 $barn_owl_mo $zfinch_proteins > barnowlMO_zfinch.paf
miniprot -t16 $barn_owl_fa $zfinch_proteins > barnowl_zfinch.paf
```

Parsing of the feature table, combination with miniprot results in PAF format and initial result exploration was conducted in [feature_paf_parsing.ipynb](feature_paf_parsing.ipynb). Alignments were filtered based on the following criteria:

- Percent identity (`p_id`): the number of matching nucleotides over the aligned nucleotides. Threshold: `> 70%`
- Percent protein covered (`p_cov`): the length of the aligned reference protein over the protein's length.  Threshold: `> 90%`

## Visualization with Sankey

Sankey visualization was conducted with R: [Sankey_protein_synteny_SupMat.Rmd](Sankey_protein_synteny_SupMat.Rmd).

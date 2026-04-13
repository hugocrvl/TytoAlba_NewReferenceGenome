# MicroFinder

Mathers *et al.* published "MicroFinder: Conserved gene-set mapping and assembly ordering for manual curation of bird microchromosomes" (2025). We used this tool to recover unplaced contigs potentially forming part of microchromosomes  to integrate them in the haplotype assemblies.In their article, they state: "*During testing we found that macrochromosome scaffolds can sometimes contain a low number of MicroFinder hits, most likely due to the presence of divergent paralogs or mis-mapping. We therefore recommend using a 5 Mb maximum scaffold size cutoff for assembly sorting.*"

We exectuted Microfinder using the container linked in the [GitHub repository](https://github.com/sanger-tol/MicroFinder).

```bash
singularity run -B `pwd`:/data docker://ghcr.io/sanger-tol/microfinder:main /data/unplaced_Fa.fa 5000000
singularity run -B `pwd`:/data docker://ghcr.io/sanger-tol/microfinder:main /data/unplaced_Mo.fa 5000000
```

MicroFinder outputs 3 files per run:

- `gff`: MiniProt alignment result (unfiltered).
- `.order.tsv`: counts file of MicroFinder hits per scaffold (filtered with threshold 70% identity).
- `ordered.fa`: MicroFinder-ordered assembly FASTA (for Hi-C subsequent analysis).

This list of 'rescued' contigs was used to further guide the assembly.

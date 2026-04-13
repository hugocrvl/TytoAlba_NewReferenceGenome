
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Exon length

get.proportion.exon.size = function(annotation.path, genome.index.path) {
  
  annotation = read.table(annotation.path, sep = "\t", fill = TRUE,  comment.char = "#", quote = "")
  colnames(annotation) = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  exons = annotation[annotation$type == "exon", ]
  genome.index = read.table(genome.index.path, header = F)
  
  proportion.exons.chr = c()
  
  for (chr in unique(exons$seqid)) {
    
    exons.chr = exons[exons$seqid == chr, ]
    exons.chr.length = exons.chr$end - exons.chr$start
    
    chr.size = genome.index$V2[genome.index$V1 == chr]
    
    proportion.exons.chr = c(proportion.exons.chr, sum(exons.chr.length) / chr.size)
    
  }
  
  output = data.frame(chromosome.size = genome.index$V2, relative.exonic.length = proportion.exons.chr)
  
  return(output)
  
}

Father.exons = get.proportion.exon.size(annotation.path= "../../13_GenomeAnnotation/3_FatherAnnotation/complete.genomic.gff",
                                               genome.index.path = "../../0_FinalFiles/Final_FatherAssembly.fa.fai")

Mother.exons = get.proportion.exon.size(annotation.path = "../../13_GenomeAnnotation/5_MotherAnnotation/complete.genomic.gff",
                                               genome.index.path = "../../0_FinalFiles/Final_MotherAssembly.fa.fai")

png("1_ExonProportion.png", res = 300, width = 4500, height = 3000)

par(mar = c(7.1, 7.1, 1.1, 1.1), mgp = c(3, 2, 0))
# Father
plot(Father.exons$chromosome.size, Father.exons$relative.exonic.length, pch = 20, col = "#156079", cex = 2, cex.axis = 2,
     xlab = "", ylab = "", xaxt = "n")
abline(lm(Father.exons$relative.exonic.length ~ Father.exons$chromosome.size), lwd = 2, col = "#156079")
# Mother
points(Mother.exons$chromosome.size, Mother.exons$relative.exonic.length, pch = 20, col = "#885A5A", cex = 2)
abline(lm(Mother.exons$relative.exonic.length ~ Mother.exons$chromosome.size), lwd = 2, col = "#885A5A")
# Axis and legends
axis(1, seq(0, max(c(Father.exons$chromosome.size, Mother.exons$chromosome.size)), 10e6), 
     labels = seq(0, max(c(Father.exons$chromosome.size, Mother.exons$chromosome.size)) / 1e6, 10), cex.axis = 2)
mtext("Chromosome size (Mb)", side = 1, line = 5, cex = 2)
mtext("Relative proportion of exons", side = 2, line = 5, cex = 2)
legend("topright", c("Paternal chromosome", "Maternal chromosome"), col = c("#156079", "#885A5A"), pch = 20, cex = 2, 
       bty = "n", inset = 0.05)

dev.off()

# Repeats

get.proportion.repeats = function(repeats.path, genome.index.path) {
  
  repeats = read.table(repeats.path, header = F)
  colnames(repeats) = c("chromosome", "start", "end")
  genome.index = read.table(genome.index.path, header = F)
  
  proportion.repeats.chr = c()
  
  chrs = unique(repeats$chromosome)
  chrs.num = as.numeric(sub("Fa_chr", "", chrs))
  chrs.sorted = chrs[order(chrs.num, na.last = T)]
  
  for (chr in chrs.sorted) {
    repeats.chr = repeats[repeats$chromosome == chr, ]
    repeats.chr.length = sum(repeats.chr$end - repeats.chr$start)
    
    genome.index.chr = genome.index$V2[genome.index$V1 == chr]
    
    proportion.repeats.chr = c(proportion.repeats.chr, repeats.chr.length / genome.index.chr)
    
  }
  
  cbind(chrs.sorted, proportion.repeats.chr)
  
  output = data.frame(chromosome.size = genome.index$V2, relative.repeats.length = proportion.repeats.chr)
  return(output)
  
}

Father.repeats = get.proportion.repeats(repeats.path = "Repeats_Sorted_Merged_Fa.txt",
                                        genome.index.path = "../../0_FinalFiles/Final_FatherAssembly.fa.fai")
Mother.repeats = get.proportion.repeats(repeats.path = "Repeats_Sorted_Merged_Mo.txt",
                                        genome.index.path = "../../0_FinalFiles/Final_MotherAssembly.fa.fai")

png("2_RepeatsProportion.png", res = 300, width = 4500, height = 3000)

par(mar = c(7.1, 7.1, 1.1, 1.1), mgp = c(3, 2, 0))
# Father
plot(Father.repeats$chromosome.size, Father.repeats$relative.repeats.length, pch = 20, col = "#156079", cex = 2, cex.axis = 2,
     xlab = "", ylab = "", xaxt = "n")
abline(lm(Father.repeats$relative.repeats.length ~ Father.repeats$chromosome.size), lwd = 2, col = "#156079")
# Mother
points(Mother.repeats$chromosome.size, Mother.repeats$relative.repeats.length, pch = 20, col = "#885A5A", cex = 2)
abline(lm(Mother.repeats$relative.repeats.length ~ Mother.repeats$chromosome.size), lwd = 2, col = "#885A5A")
# Axis and legends
axis(1, seq(0, max(c(Father.repeats$chromosome.size, Mother.repeats$chromosome.size)), 10e6), 
     labels = seq(0, max(c(Father.repeats$chromosome.size, Mother.repeats$chromosome.size)) / 1e6, 10), cex.axis = 2)
mtext("Chromosome size (Mb)", side = 1, line = 5, cex = 2)
mtext("Relative proportion of repeats", side = 2, line = 5, cex = 2)
legend("topright", c("Paternal chromosome", "Maternal chromosome"), col = c("#156079", "#885A5A"), pch = 20, cex = 2, 
       bty = "n", inset = 0.05)

dev.off()

# GC

get.GC = function(GC.path, genome.index.path) {
  
  GC = read.table(GC.path, header = F)
  genome.index = read.table(genome.index.path, header = F)
  
  output = data.frame(chromosome.size = genome.index$V2, GC = GC$V1)
  return(output)
  
}

GC.Father = get.GC(GC.path = "GC_Father.txt",
                   genome.index.path = "../../0_FinalFiles/Final_FatherAssembly.fa.fai")
GC.Mother = get.GC(GC.path = "GC_Mother.txt",
                   genome.index.path = "../../0_FinalFiles/Final_MotherAssembly.fa.fai")

png("3_GCProportion.png", res = 300, width = 4500, height = 3000)

par(mar = c(7.1, 7.1, 1.1, 1.1), mgp = c(3, 2, 0))
# Father
plot(GC.Father$chromosome.size, GC.Father$GC, pch = 20, col = "#156079", cex = 2, cex.axis = 2,
     xlab = "", ylab = "", xaxt = "n")
abline(lm(GC.Father$GC ~ GC.Father$chromosome.size), lwd = 2, col = "#156079")
# Mother
points(GC.Mother$chromosome.size, GC.Mother$GC, pch = 20, col = "#885A5A", cex = 2)
abline(lm(GC.Mother$GC ~ GC.Mother$chromosome.size), lwd = 2, col = "#885A5A")
# Axis and legends
axis(1, seq(0, max(c(GC.Father$chromosome.size, GC.Mother$chromosome.size)), 10e6), 
     labels = seq(0, max(c(GC.Father$chromosome.size, GC.Mother$chromosome.size)) / 1e6, 10), cex.axis = 2)
mtext("Chromosome size (Mb)", side = 1, line = 5, cex = 2)
mtext("Relative proportion of GC", side = 2, line = 5, cex = 2)
legend("topright", c("Paternal chromosome", "Maternal chromosome"), col = c("#156079", "#885A5A"), pch = 20, cex = 2, 
       bty = "n", inset = 0.05)

dev.off()

# New genes

library(dplyr)

get.best.CDS = function(CDS.alignment.path) {
  
  alignment = read.table(CDS.alignment.path, header = F)
  colnames(alignment) = c("QueryName", "QueryLength", "QueryStart", "QueryEnd", "TargetName", "MappingQuality")
  length.alignment = alignment$QueryEnd - alignment$QueryStart
  # Keeping fully aligned genes, with MAPQ = 60
  filtering = alignment$QueryLength == length.alignment & alignment$MappingQuality == 60
  
  common = data.frame(chr = sub(".*_chr([A-Za-z0-9]+).*", "\\1", alignment$QueryName[filtering]), alignment[filtering, ], FILTER = "Common")
  new = data.frame(chr = sub(".*_chr([A-Za-z0-9]+).*", "\\1", alignment$QueryName[!filtering]), alignment[!filtering, ], FILTER = "New")
  
  new_best = new %>%
    group_by(QueryName) %>%
    slice_max(order_by = QueryLength, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  output = rbind(common, new_best)
  
  return(output)
  
}

get.CDS.figure = function(CDS.alignment.path, png.name, colour.common, colour.new) {
  
  CDS.alignment = get.best.CDS(CDS.alignment.path)

  # Absolute number
  
  total.count = table(CDS.alignment$chr)
  total.count = total.count[order(as.numeric(names(total.count)), na.last = T)]
  chr.name = names(total.count)
  
  tmp = table(CDS.alignment$chr[CDS.alignment$FILTER == "New"])
  tmp = tmp[order(as.numeric(names(tmp)), na.last = T)]
  new.count = rep(0, 46)
  new.count[-(as.numeric(setdiff(chr.name, names(tmp))))] = tmp
  names(new.count) = chr.name

  # Proportion
  
  proportion.new = new.count / total.count[names(new.count)]

  # Reorder the vectors
  
  # chr.levels = c(seq(1, 46)[seq(1, 46) != 40], "Z")
  # proportion.new.reorder = proportion.new[chr.levels]
  # new.count.reorder = new.count[chr.levels]
  # total.count.reorder = total.count[chr.levels]
  
  png(png.name, res = 300, width = 4500, height = 4250)
  
  layout(matrix(c(1, 2), ncol = 1),
         heights = c(0.35, 0.65), widths = c(1))
  
  par(mar = c(5.1, 9.1, 3.1, 1.1), mgp = c(3, 2, 0))
  
  barplot(rep(1, 46), width = 1, yaxt = "n", col = colour.common, cex.axis = 2)
  barplot(proportion.new, col = colour.new, add = T, las = 2, cex.axis = 2)
  mtext("Relative proportion of CDS", side = 2, line = 7, cex = 2)
  
  par(mar = c(7.1, 9.1, 1.1, 1.1), mgp = c(3, 2, 0))
  
  barplot(total.count, col = colour.common, las = 2, cex.axis = 2, xaxt = "n")
  barplot(new.count, col = colour.new, las = 2, yaxt = "n",
          add = T)
  mtext("Absolute number of CDS", side = 2, line = 7, cex = 2)
  mtext("Chromosomes", side = 1, line = 5, cex = 2)
  
  dev.off()
  
}

get.CDS.figure(CDS.alignment.path = "../../13_GenomeAnnotation/7.2_AlignedCDS2ASSEMBLY/FilteredFatherCDS_Aligned2UNZIP_Subset.txt",
               png.name = "4.1_GenesFather.png",
               colour.common = "#156079",
               colour.new = "#95b1af")
               
get.CDS.figure(CDS.alignment.path = "../../13_GenomeAnnotation/7.2_AlignedCDS2ASSEMBLY/FilteredMotherCDS_Aligned2UNZIP_Subset.txt",
               png.name = "4.2_GenesMother.png",
               colour.common = "#885A5A",
               colour.new = "#fed9b7")

# Getting the names of the new CDS only

Father = get.best.CDS(
  CDS.alignment.path = "../../13_GenomeAnnotation/7.2_AlignedCDS2ASSEMBLY/FilteredFatherCDS_Aligned2UNZIP_Subset.txt"
  )
write.table(Father$QueryName[Father$FILTER == "New"], "../../13_GenomeAnnotation/6_ProcessingOutput/Father_NameCDS.txt", quote = F, col.names = F, row.names = F)

Mother = get.best.CDS(
  CDS.alignment.path = "../../13_GenomeAnnotation/7.2_AlignedCDS2ASSEMBLY/FilteredMotherCDS_Aligned2UNZIP_Subset.txt"
)
write.table(Mother$QueryName[Mother$FILTER == "New"], "../../13_GenomeAnnotation/6_ProcessingOutput/Mother_NameCDS.txt", quote = F, col.names = F, row.names = F)

################

# Creation of the data.frame for the chromosome categorisation
# I previously computed all the *.exons, *.repeats and GC.* data.frame

index.father = read.table("../../0_FinalFiles/Final_FatherAssembly.fa.fai", header = F)
chromosomes.Father = data.frame(chr = index.father$V1,
                                size = Father.exons$chromosome.size, 
                                exons = Father.exons$relative.exonic.length,
                                repeats = Father.repeats$relative.repeats.length,
                                GC = GC.Father$GC)
write.table(chromosomes.Father, "../../15_ChromosomeCategories/ProportionsFather.txt", quote = F, row.names = F, col.names = T)

index.mother = read.table("../../0_FinalFiles/Final_MotherAssembly.fa.fai", header = F)
chromosomes.Mother = data.frame(chr = index.mother$V1,
                                size = Mother.exons$chromosome.size, 
                                exons = Mother.exons$relative.exonic.length,
                                repeats = Mother.repeats$relative.repeats.length,
                                GC = GC.Mother$GC)
write.table(chromosomes.Mother, "../../15_ChromosomeCategories/ProportionsMother.txt", quote = F, row.names = F, col.names = T)

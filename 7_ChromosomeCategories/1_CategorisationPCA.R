
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Function to get the final data.frame

get.zscores = function(df) {
  
  df.zscores = data.frame(chr = df$chr)
  for (i in 2:ncol(df)) {
    
    col = df[, i]
    z.scores = (col - mean(col)) / sd(col)
    df.zscores = cbind(df.zscores, z.scores)
    colnames(df.zscores)[i] = colnames(df)[i]
    
  }
  
  return(df.zscores)
  
}

get.data.frame.PCA = function(proportions.path, genes.path, repeats.path) {
  
  proportions = read.table(proportions.path, header = T)
  genes = read.table(genes.path, header = T)
  genes = genes[order(as.numeric(sub(".*_chr(.*)", "\\1", genes$Chromosome)), na.last = T), ]
  genes$Genes = genes$Genes / proportions$size
  genes$Exons_unique = genes$Genes / proportions$size
  repeats = read.table(repeats.path, header = T)
  repeats = repeats[order(as.numeric(sub(".*_chr(.*)", "\\1", repeats$Chromosome)), na.last = T), ]
  repeats.relative = repeats[, -1] / proportions$size
  
  characteristics = cbind(proportions[, -4], genes[, -1], repeats.relative)
  characteristics.zscores = get.zscores(characteristics)
  
  return(characteristics.zscores)
  
}

# Father

characteristics.Father.zscores = get.data.frame.PCA(
  proportions.path = "ProportionsFather.txt",
  genes.path = "GenesFeaturesFather.txt",
  repeats.path = "RepeatsTypesFather.txt"
  )

PCA.Father = prcomp(characteristics.Father.zscores[, -c(1)])
var.explained.Father = PCA.Father$sdev^2 / sum(PCA.Father$sdev^2)

# Mother

characteristics.Mother.zscores = get.data.frame.PCA(
  proportions.path = "ProportionsMother.txt",
  genes.path = "GenesFeaturesMother.txt",
  repeats.path = "RepeatsTypesMother.txt"
)

PCA.Mother = prcomp(characteristics.Mother.zscores[, -c(1)])
var.explained.Mother = PCA.Mother$sdev^2 / sum(PCA.Mother$sdev^2)

# Making figure

new.chrs.Fa = match(c("Fa_chr33", "Fa_chr40", "Fa_chr42", "Fa_chr43", "Fa_chr44", "Fa_chr45"), 
                 characteristics.Father.zscores$chr)
new.chrs.Mo = match(c("Mo_chr33", "Mo_chr40", "Mo_chr42", "Mo_chr43", "Mo_chr44", "Mo_chr45"), 
                    characteristics.Mother.zscores$chr)

# Saving the non-z-scores transformed versions
# wb = openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb, "PaternalHaplome")
# openxlsx::addWorksheet(wb, "MaternalHaplome")
# 
# openxlsx::writeData(wb, sheet = "PaternalHaplome", characteristics.Father.zscores)
# openxlsx::writeData(wb, sheet = "MaternalHaplome", characteristics.Mother.zscores)
# 
# openxlsx::setColWidths(wb, sheet = "PaternalHaplome", cols = 1:ncol(characteristics.Father.zscores), widths = "auto")
# openxlsx::setColWidths(wb, sheet = "MaternalHaplome", cols = 1:ncol(characteristics.Mother.zscores), widths = "auto")
# 
# openxlsx::saveWorkbook(wb, "SupplementaryTable_ChromosomeFeatures.xlsx", overwrite = TRUE)

png("1.1_PCA_ChromosomesFeatures.png", height = 2250, width = 4500, res = 300)

layout(matrix(c(1, 2), ncol = 2))

par(mar = c(7.1, 7.1, 1.1, 1.1))
par(mgp = c(3, 2, 0))

plot(PCA.Father$x[, 1], PCA.Father$x[, 2], pch = 20, col = "#156079", 
     xlab = "", ylab = "", cex = 2, cex.axis = 2)
mtext(paste0("PC1 ", round(var.explained.Father[1]*100, 2), "%"), side = 1, line = 5, cex = 2)
mtext(paste0("PC2 ", round(var.explained.Father[2]*100, 2), "%"), side = 2, line = 5, cex = 2)
points(PCA.Father$x[new.chrs.Fa, 1], PCA.Father$x[new.chrs.Fa, 2], pch = 20, 
       col = "#95b1af", cex = 2)
text(PCA.Father$x[, 1], PCA.Father$x[, 2], gsub("Fa_chr", "", characteristics.Father.zscores$chr),
     pos = 3, cex = 1)

par(mar = c(7.1, 7.1, 1.1, 1.1))
par(mgp = c(3, 2, 0))

plot(PCA.Mother$x[, 1], PCA.Mother$x[, 2], pch = 20, col = "#885a5a", 
     xlab = "", ylab = "", cex = 2, cex.axis = 2)
mtext(paste0("PC1 ", round(var.explained.Mother[1]*100, 2), "%"), side = 1, line = 5, cex = 2)
mtext(paste0("PC2 ", round(var.explained.Mother[2]*100, 2), "%"), side = 2, line = 5, cex = 2)
points(PCA.Mother$x[new.chrs.Mo, 1], PCA.Mother$x[new.chrs.Mo, 2], pch = 20, 
       col = "#fed9b7", cex = 2)
text(PCA.Mother$x[, 1], PCA.Mother$x[, 2], gsub("Mo_chr", "", characteristics.Mother.zscores$chr),
     pos = 3, cex = 1)

dev.off()

# Overlap

plot(PCA.Father$x[, 1], PCA.Father$x[, 2], pch = 20, col = "#156079", 
     xlab = "PC1", ylab = "PC2", cex = 2)
text(PCA.Father$x[, 1], PCA.Father$x[, 2], gsub("Fa_chr", "", characteristics.Father.zscores$chr),
     pos = 4, cex = 0.7, col = "#156079")

points(PCA.Mother$x[, 1], PCA.Mother$x[, 2], pch = 20, col = "#885a5a",
       cex = 2)
text(PCA.Mother$x[, 1], PCA.Mother$x[, 2], gsub("Mo_chr", "", characteristics.Mother.zscores$chr),
     pos = 4, cex = 0.7, col = "#885a5a")

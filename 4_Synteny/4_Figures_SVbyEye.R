
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(SVbyEye)

dir.create("4.1_Alignment_Figures", showWarnings = F)

for (path in list.files("2.2_Aligned_LGs/", full.names = T)) {
  
  print(path)
  paf.table = readPaf(paf.file = path, include.paf.tags = TRUE, restrict.paf.tags = "cg")
  
  tmp = gsub("Aligned_", "", basename(path))
  LG = gsub(".paf", "", tmp)
  
  png(paste0("4.1_Alignment_Figures/Aligned_", LG, ".png"), res = 300, height = 3000, width = 4500)
  figure = plotMiro(paf.table = paf.table, color.by = "direction", color.palette = c("-" = "#696863", "+" = "#C0BDA5")) + 
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "#FFFFFF"),
      plot.background = ggplot2::element_rect(fill = "#FFFFFF")
    ) 
  plot(figure)
  dev.off()
}

#####################
#####################
#####################

library("pafr")
paf = read_paf("3.2_AlignmentWholeGenome.paf")
order.query = sort(unique(paf$qname))
order.target = sort(unique(paf$tname))
orders = list(order.query, order.target)

png("4.1_Alignment_Figures/AlignedWholeGenome.png", res = 300, height = 2500, width = 2500)

dotplot(paf, order_by = "provided", ordering = orders, dashes = T,
        alignment_colour = "#0B4F6C", xlab = "Mother LGs", ylab = "Father LGs",
        line_size = 2) + theme_bw() + theme(panel.background = element_rect(fill = "#D9D8D3"),
                                          plot.background = element_rect(fill = "#D9D8D3"))

dev.off()

#####################
#####################
#####################


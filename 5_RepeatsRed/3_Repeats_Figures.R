
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

####################################################

get.middle.position = function(start.vector, end.vector) {
  
  mid = start.vector + ((end.vector - start.vector) / 2)
  return(mid)
}

get.absolute.position = function(relative.position.vector, scaffold.vector) {
  
  scaffold.vector = factor(scaffold.vector, levels = unique(scaffold.vector))
  # posmin = tapply(position.vector, scaffold.vector, min)
  posmax = tapply(relative.position.vector, scaffold.vector, max)
  posshift = head(c(0, cumsum(posmax)), -1) 
  names(posshift) = unique(scaffold.vector)
  abs_pos = c()
  
  for (i in 1:length(relative.position.vector)) {
    
    chr = scaffold.vector[i]
    pos = relative.position.vector[i]
    abs_pos[i] = posshift[as.character(chr)] + pos
    
  }
  
  return(abs_pos)
  
}

get.thick.mark = function(relative.position.vector, scaffold.vector) {
  scaffold.vector = factor(scaffold.vector, levels = unique(scaffold.vector))
  posmin = tapply(relative.position.vector, scaffold.vector, min);
  posmax = tapply(relative.position.vector, scaffold.vector, max);
  posshift = head(c(0, cumsum(posmax)), -1); 
  names(posshift) = unique(scaffold.vector);
  mid = posmin + ((posmax-posmin)/2) 
  mid = mid + posshift
  return(mid)
}

add.alpha = function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

get.scaffold.colour = function(scaffold.vector) {
  scaffold.factor = factor(scaffold.vector, levels = unique(scaffold.vector))
  if (length(levels(scaffold.factor)) %% 2 == 0) {
    colour = rep(c("#999999", "#EAEAEA"), length(levels(scaffold.factor)) %/% 2)[scaffold.factor]
  } else {
    colour.name = c(rep(c("#999999", "#EAEAEA"), length(levels(scaffold.factor)) %/% 2), "#999999")
    colour = rep(colour.name)[scaffold.factor]
  }
  return(colour)
}

####################################################

# Father haplotype

repeats.Fa = read.table("all_Fa.100000bp.pct.bedgraph", header = F)
colnames(repeats.Fa) = c("chrom", "start", "end", "pct_repeats")

mid.Fa = get.middle.position(repeats.Fa$start, repeats.Fa$end)
abs.mid.Fa = get.absolute.position(mid.Fa, repeats.Fa$chrom)
thick.mark.Fa = get.thick.mark(mid.Fa, repeats.Fa$chrom)
colours.Fa = get.scaffold.colour(repeats.Fa$chrom)

png("3.1_Repeats_Fa.png", res = 300, height = 3000, width = 6000)
par(mgp = c(4,2,0))
par(mar = c(7.1, 10.1, 3.1, 1.1))
plot(abs.mid.Fa, repeats.Fa$pct_repeats, pch = 20, col = colours.Fa, xaxt = "n",
     xlab = "", ylab = "", cex.lab = 2, cex.axis = 2, cex.main = 1.5, cex = 2,
     main = paste0("Father haplotype, whole genome average = ", round(mean(repeats.Fa$pct_repeats), 2)))
axis(1, at = thick.mark.Fa, labels = names(thick.mark.Fa), las = 2, cex.axis = 1.25)
mtext("Percentage of repeats", side = 2, line = 5, cex = 2)
dev.off()


# Mother haplotype

repeats.Mo = read.table("all_Mo.100000bp.pct.bedgraph", header = F)
colnames(repeats.Mo) = c("chrom", "start", "end", "pct_repeats")

mid.Mo = get.middle.position(repeats.Mo$start, repeats.Mo$end)
abs.mid.Mo = get.absolute.position(mid.Mo, repeats.Mo$chrom)
thick.mark.Mo = get.thick.mark(mid.Mo, repeats.Mo$chrom)
colours.Mo = get.scaffold.colour(repeats.Mo$chrom)

png("3.2_Repeats_Mo.png", res = 300, height = 3000, width = 6000)
par(mgp = c(4,2,0))
par(mar = c(7.1, 10.1, 3.1, 1.1))
plot(abs.mid.Mo, repeats.Mo$pct_repeats, pch = 20, col = colours.Fa, xaxt = "n",
     xlab = "", ylab = "", cex.lab = 2, cex.axis = 2, cex.main = 1.5, cex = 2,
     main = paste0("Mother haplotype, whole genome average = ", round(mean(repeats.Mo$pct_repeats), 2)))
axis(1, at = thick.mark.Mo, labels = names(thick.mark.Mo), las = 2, cex.axis = 1.25)
mtext("Percentage of repeats", side = 2, line = 5, cex = 2)
dev.off()
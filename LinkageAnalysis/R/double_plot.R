# this function plots the heatmap of p values
heatmap_pval <- function(sig, type, test, cutoff_pair, mask, fns) {
  # break points for color scheme of heatmap
  rg <- range(-log10(sig[[type]]))
  if (rg[2] == Inf) {
    rg[2] <- -log10(cutoff_pair)
  }
  mid <- round(min(-log10(cutoff_pair), rg[2]/2))
  
  # the significant pairs should have better resolution
  a <- seq(from = 0, to = mid, by = mid/2)
  b <- seq(from = mid, to = max(ceiling(rg[2]), -log10(cutoff_pair)), length.out = 10)
  c <- unique(c(a, b))
  c <- c[order(c)]
  
  if (all(sig[[type]] == 1)) {
    report("e", "All p values are equal to 1! \n    Maybe due to too few number of mice!", 
      fns$log_file)
  }
  
  heatmap.2(-log10(sig[[type]][mask, mask]), dendrogram = "none", trace = "none", 
    Rowv = F, Colv = F, cexRow = 0.7, cexCol = 0.7, main = paste("Heatmap of P values\n (", 
      type, test, ")", sep = " "), lmat = rbind(c(2, 3, 4), c(0, 1, 1)), lhei = c(1, 
      5), lwid = c(0.5, 3, 1), breaks = c, keysize = 0.5, density.info = "none")
}

# this function plots the distribution plot of genotype vs. phenotype
double_distrib <- function(flag, data, input, type, pair) {
  title <- paste(pair[1, 1], pair[1, 2], "\n", pair[2, 1], pair[2, 2], "\n", flag)
  x_ticks <- c("0 0", "0 1", "0 2", "1 0", "1 1", "1 2", "2 0", "2 1", "2 2")
  x_labs <- c("REF\nREF", "REF\nHET", "REF\nVAR", "HET\nREF", "HET\nHET", "HET\nVAR", 
    "VAR\nREF", "VAR\nHET", "VAR\nVAR")
  
  G3_type <- paste(data$gt1, data$gt2)  # mark the # G3 mice with each genotype combination
  for (i in 1:9) {
    x_labs[i] <- paste(x_labs[i], "\n(", sum(G3_type == x_ticks[i]), ")", sep = "")
  }
  
  if (input$bin == FALSE) {
    # continuous mode plot distribution and x axis
    x_noise <- runif(length(data$gt1), 0, 0.4) - 0.2
    plot(data$gt1 * 3 + data$gt2 + x_noise, data$pt, xlab = "", ylab = "", xaxt = "n", 
      yaxt = "n", main = title, xlim = c(0, 9.5), pch = 19, col = as.numeric(data$mother), 
      cex.main = 0.9)
    axis(side = 1, at = c(0:8), labels = x_labs, tick = FALSE, cex.axis = 0.6)
    axis(side = 2, las = 1)
    
    # legend
    rg <- range(data$pt)
    tmp <- seq(from = rg[1], to = rg[2], length.out = length(unique(data$mother)) + 
      2)
    text(x = 9, y = tmp[-c(1, length(tmp))], labels = unique(data$mother), col = as.numeric(unique(data$mother)))
  } else {
    # binary mode
    
    plot(1:10, 1:10, xlim = c(0, 1), ylim = c(0, 1), xaxt = "n", yaxt = "n", 
      main = title, axes = F, xlab = "", ylab = "Affected G3 mice%", type = "n", 
      cex.main = 0.9)
    axis(side = 2, at = c(0:5)/5, labels = c(0:5) * 20, las = 1)  # y axis (left ticks)
    left <- 0
    
    for (i in 1:length(x_ticks)) {
      tot <- sum(paste(data$gt1, data$gt2) == x_ticks[i])
      if (tot == 0) {
        next
      }
      aff <- sum(paste(data$gt1, data$gt2) == x_ticks[i] & data$pt == "AFFECTED")
      
      rect(xleft = left + tot/dim(data)[1] * 0.1, xright = left + tot/dim(data)[1] * 
        0.9, ybottom = 0, ytop = 1, col = "gray91")  # total count
      rect(xleft = left + tot/dim(data)[1] * 0.1, xright = left + tot/dim(data)[1] * 
        0.9, ybottom = 0, ytop = aff/tot, col = "gray21")  # mutant count
      axis(side = 1, at = left + tot/dim(data)[1] * 0.5, labels = x_labs[i], 
        cex.axis = 0.5, tick = FALSE)
      left <- left + tot/dim(data)[1]
    }
  }
} 

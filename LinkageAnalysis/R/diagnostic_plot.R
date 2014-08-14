# This function only works if the phenotype scores are continuous 1) plots a
# density plot with histogram 2) plots a CDF plots of phenotype score
diagnostic <- function(phenotype, bin, log_file, unconverted) {
  # no plotting if the phenotype is binary and not converted from a continuous
  # version
  if (bin == TRUE && is.null(unconverted)) {
    return(0)
  }
  par(mar = c(5, 5, 5, 5))

  # in this function, phenotype is changed. The continuous scores in unconverted
  # may be used to cover the binary version in phenotype cutoff is defined only
  # when phenotype is converted
  if (!is.null(unconverted)) {
    phenotype$phenotype <- unconverted$phenotype
    par(mfrow = c(1, 2))
    cutoff <- c(unconverted$break_point, 1e+08)
  } else {
    par(mfrow = c(1, 1))
    cutoff <- NULL
  }

  # plot density plot with histogram
  rg <- range(phenotype$phenotype)  # set bin length of histogram
  bin_step <- (rg[2] - rg[1])/20
  rg[1] <- floor(rg[1]/bin_step) * bin_step
  rg[2] <- ceiling(rg[2]/bin_step) * bin_step

  hist(phenotype$phenotype, breaks = seq(rg[1], rg[2], by = bin_step), xlab = "Phenotype scores",
    main = "Histogram of phenotype scores", cex.lab = 2, cex.main = 2, border = "black",
    col = "red", prob = T, cex.axis = 2)
  den <- density(phenotype$phenotype, from = rg[1], to = rg[2])
  lines(den, lwd = 4, col = "gold1")
  if (!is.null(cutoff)) {
    for (i in 1:length(cutoff)) {
      lines(x = c(cutoff[i], cutoff[i]), y = c(0, max(den$y)), lwd = 4, col = "blue")
    }
  }

  # plot CDF plots
  if (!is.null(cutoff)) {
    pheno_scores <- phenotype$phenotype
    pheno_scores <- pheno_scores[order(pheno_scores)]  # in case the phenotype scores are not sorted

    pheno_col <- rep(0, length(pheno_scores))  # assign colors
    for (i in 1:length(cutoff)) {
      pheno_col[pheno_scores < cutoff[i]] <- pheno_col[pheno_scores < cutoff[i]] +
        1
    }

    plot(1:length(pheno_scores), pheno_scores, pch = 19, xlab = "G3 mouse", main = "Phenotype score CDF plot",
      ylab = "Phenotype scores", cex.lab = 2, cex.axis = 2, col = pheno_col,
      cex.main = 2, ylim = c(rg[1], rg[2]))

    # draw legend
    for (i in 1:length(cutoff)) {
      x_pos <- 5/6 * length(pheno_scores)
      y_pos <- i/10 * (rg[2] - rg[1]) + rg[1]
      text(x = x_pos, y = y_pos, lab = paste("Group", i), col = length(cutoff) +
        1 - i, font = 2, cex = 2)
    }

    par(mfrow = c(1, 1))
  }
}


plot.hist <- function(x, main = "") {
  x <- data.frame(x = as.numeric(x))
  g <- ggplot(x, aes(x = x)) +
      geom_histogram(aes(y = ..density.., fill = ..count..)) +
      geom_density(col = "gold1", size = 1.5) + xlab("Phenotype scores") + ylab("Density") + ggtitle(main)
  print(g)
  g
}

# this function draws scattorplots of phenotypes scores (continuous variable) or
# table plots of phenotypes (binary variable).
distrib_plot <- function(data, genes_i, bin) {
  if (bin == TRUE) {
    plot(factor(data$gt), data$pt, xlab = "", ylab = "Affected G3 mice%", axes = FALSE,
         main = paste(genes_i[, c("Gene", "Coordination")], collapse = "\n"))  # table plot

    exist_gt <- as.numeric(as.vector(unique(factor(data$gt))))
    exist_gt <- exist_gt[order(exist_gt)]
    exist_gt <- as.vector(factor(exist_gt, levels = c(0, 1, 2), labels = c("REF",
                                                                    "HET", "VAR")))
    if (length(exist_gt) == 3) {
      axis(side = 1, at = c(0.2, 0.5, 0.8), labels = c("REF", "HET", "VAR"),
           tick = F)  # x axis
    } else {
      # some genes don't have all three genotypes
      axis(side = 1, at = c(0.25, 0.75), labels = exist_gt, tick = F)
    }

    axis(side = 2, at = c(0:5)/5, labels = c(0:5) * 20, las = 1)  # y axis (left ticks)
  } else {
    data$gt <- data$gt + runif(dim(data)[1], -0.25, 0.25)  # add noise to data points
    leg_pos <- seq(from = min(data$pt), to = max(data$pt), by = (max(data$pt) -
                                            min(data$pt))/(length(unique(data$mother)) + 1))[2:(length(unique(data$mother)) +
                                                                                                1)]
    legends <- data.frame(pt = leg_pos, sex = 1, gt = 2.7, mother = unique(data$mother))  # fake data frame to serve as legend
    data <- rbind(legends, data)

    plot(data$gt, data$pt, xaxt = "n", pch = 19, col = data$mother, xlim = c(-0.5,
                                                                        3.5), xlab = "", ylab = "", main = paste(genes_i[, c("Gene", "Coordination")],
                                                                                                        collapse = "\n"))
    axis(side = 1, at = c(0, 1, 2), labels = c("REF", "HET", "VAR"))
    text(x = 2.7, y = legends$pt, labels = legends$mother, pos = 4)  # label legends
  }
}

plot.distribution <- function(pheno, pheno.name, main = "") {
  stopifnot(!is.null(pheno))
  stopifnot(!is.null(pheno$gt))
  stopifnot(!all(is.na(pheno$gt)))
  stopifnot(!is.null(pheno[, pheno.name]))
  stopifnot(!all(is.na(pheno[,pheno.name])))
  gt <- mother <- NULL ## bypass CRAN check
  # library(ggplot2)
  pheno$family = pheno$fid
  pheno$y <- pheno[, pheno.name]
  pheno$y <- as.numeric(pheno$y)
  ##pheno$pch <- rep_len(x = seq(6), length.out = nrow(pheno))
  y <- NULL ## bypass CRAN check
  g <- ggplot(pheno, aes(x = gt, y = y, col = family)) +
    geom_point(position = position_jitter(width = 0.2)) +
      xlab("Genotype") + ylab("Phenotype") +
        xlim(-0.5, 2.5) + ggtitle(main)
  return(g)
}

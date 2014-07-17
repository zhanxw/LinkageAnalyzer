# calculate penetrance and semidominance
get_genetics <- function(genes, phenotype, genotype, bin, unconverted) {
  genetics <- matrix(data = NaN, ncol = 4, dim(genes)[1])
  colnames(genetics) <- c("Penetrance_REF", "Penetrance_HET", "Penetrance_VAR",
    "Semidominance")

  pt <- phenotype$phenotype  # pt is the original phenotype scores (if converted)
  if (!is.null(unconverted)) {
    pt <- unconverted$phenotype
  }

  # calculate penetrance. continuous phenotype will all lead to NaN calculate
  # semidominance. skip binary phenotype scores if continuous version is not
  # available
  for (i in 1:dim(genes)[1]) {
    if (bin == TRUE) {
      # penetrance
      genetics[i, "Penetrance_REF"] <-
        sum(genotype[i, ] == "REF" & phenotype$phenotype == "AFFECTED")/sum(genotype[i, ] == "REF")
      genetics[i, "Penetrance_HET"] <-
        sum(genotype[i, ] == "HET" & phenotype$phenotype == "AFFECTED")/sum(genotype[i, ] == "HET")
      genetics[i, "Penetrance_VAR"] <-
        sum(genotype[i, ] == "VAR" & phenotype$phenotype == "AFFECTED")/sum(genotype[i, ] == "VAR")
    }

    if (bin == TRUE && is.null(unconverted)) {
      next
    }
    tmp <- (mean(pt[genotype[i, ] == "REF"]) - mean(pt[genotype[i, ] == "HET"]))/
      (mean(pt[genotype[i,] == "REF"]) - mean(pt[genotype[i, ] == "VAR"]))
    if (!is.na(tmp)) {
      if (tmp > 1) {
        tmp <- 1
      }
      if (tmp < 0) {
        tmp <- 0
      }
    }
    genetics[i, "Semidominance"] <- tmp
  }

  genetics[is.na(genetics)] <- NA
  genetics <- round(genetics, digits = 3)
  return(genetics)
}

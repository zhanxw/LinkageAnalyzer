##  ====================================================================================================================================================
##  |  LinkageAnalysis License                                                                                                                         |
##  |  ----------------------------------------------------------------------------------------------------------------------------------------------  |
##  |  a.   Copyright ©2014, The University of Texas Southwestern Medical Center.  All rights reserved; and                                            |
##  |  b.   This software and any related documentation constitutes published and/or unpublished works and may contain valuable trade secrets and      |
##  |       proprietary information belonging to The University of Texas Southwestern Medical Center (UT SOUTHWESTERN).  None of the foregoing         |
##  |       material may be copied, duplicated or disclosed without the express written permission of UT SOUTHWESTERN.  IN NO EVENT SHALL UT           |
##  |       SOUTHWESTERN BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING   |
##  |       OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF UT SOUTHWESTERN HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.         |
##  |       UT SOUTHWESTERN SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND        |
##  |       FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". UT          |
##  |       SOUTHWESTERN HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.                                   |
##  |  c.   This software contains copyrighted materials from R-package, ggplot2, gplots, gridExtra, lme4, logging, mclust, plyr and stringr.          |
##  |       Corresponding terms and conditions apply.                                                                                                  |
##  ====================================================================================================================================================    
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

calc.genetic <- function(ret, geno, pheno, pheno.name, isBinary) {
  if (ncol(geno) == 0) {
    ret$REF <- ret$HET <- ret$VAR <- rep(0, nrow(geno))
    return(ret)
  }

  ret$REF <- rowSums(geno == 0, na.rm = TRUE)
  ret$HET <- rowSums(geno == 1, na.rm = TRUE)
  ret$VAR <- rowSums(geno == 2, na.rm = TRUE)

  if (isBinary) {
    ret$Penetrance_REF <- apply(geno, 1, function(x) {idx <- x==0; mean(pheno[idx, pheno.name] == "AFFECTED", na.rm = TRUE)})
    ret$Penetrance_HET <- apply(geno, 1, function(x) {idx <- x==1; mean(pheno[idx, pheno.name] == "AFFECTED", na.rm = TRUE)})
    ret$Penetrance_VAR <- apply(geno, 1, function(x) {idx <- x==2; mean(pheno[idx, pheno.name] == "AFFECTED", na.rm = TRUE)})
  }
  ret$Semidominance <- apply(geno, 1, function(x) {
    idx <- x == 0; m0 <- mean(as.numeric(pheno[idx, pheno.name]), na.rm = TRUE)
    idx <- x == 1; m1 <- mean(as.numeric(pheno[idx, pheno.name]), na.rm = TRUE)
    idx <- x == 2; m2 <- mean(as.numeric(pheno[idx, pheno.name]), na.rm = TRUE)
    ## cat(m0, m1,m2, "\n")
    ret <- (m0 - m1) / (m0 - m2)
    if (!is.na(ret)) {
      if (ret > 1) { ret <- 1}
      if (ret < 0) { ret <- 0}
    }
    ret
  })
  ret
}

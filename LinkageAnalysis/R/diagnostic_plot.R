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

##  ====================================================================================================================================================
##  |  This file is part of LinkageAnalysis.													       |
##  |																		       |
##  |  LinkageAnalysis is free software: you can redistribute it and/or modify									       |
##  |  it under the terms of the GNU General Public License as published by									       |
##  |  the Free Software Foundation, either version 3 of the License, or									       |
##  |  (at your option) any later version.													       |
##  |																		       |
##  |  LinkageAnalysis is distributed in the hope that it will be useful,									       |
##  |  but WITHOUT ANY WARRANTY; without even the implied warranty of										       |
##  |  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the										       |
##  |  GNU General Public License for more details.												       |
##  |																		       |
##  |  You should have received a copy of the GNU General Public License									       |
##  |  along with LinkageAnalysis.  If not, see <http://www.gnu.org/licenses/>.									       |
##  ====================================================================================================================================================


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
  ..density.. <- ..count.. <- NULL ## bypass CRAN check
  g <- ggplot(x, aes(x = x)) +
      geom_histogram(aes(y = ..density.., fill = ..count..)) +
      geom_density(col = "gold1", size = 1.5) + xlab("Phenotype scores") + ylab("Density") + ggtitle(main)
  print(g)
  g
}

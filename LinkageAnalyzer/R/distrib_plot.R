##  ====================================================================================================================================================
##  |  LinkageAnalyzer Copyright                                                                                                                       |
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
##  |  This file is part of LinkageAnalyzer.                                                                                                           |
##  |                                                                                                                                                  |
##  |  LinkageAnalyzer is free software: you can redistribute it and/or modify                                                                         |
##  |  it under the terms of the GNU General Public License as published by                                                                            |
##  |  the Free Software Foundation, either version 3 of the License, or                                                                               |
##  |  (at your option) any later version.                                                                                                             |
##  |                                                                                                                                                  |
##  |  LinkageAnalyzer is distributed in the hope that it will be useful,                                                                              |
##  |  but WITHOUT ANY WARRANTY; without even the implied warranty of                                                                                  |
##  |  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                                                   |
##  |  GNU General Public License for more details.                                                                                                    |
##  |                                                                                                                                                  |
##  |  You should have received a copy of the GNU General Public License                                                                               |
##  |  along with LinkageAnalyzer.  If not, see <http://www.gnu.org/licenses/>.                                                                        |
##  ====================================================================================================================================================


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

# arrange plots into multiple subfigures and save to a PDF file
# dist.plots list of ggplot2 objects
save.dist.plot <- function(dist.plots, nrow, ncol, fn) {
  nplots <- length(dist.plots)
  loginfo("%d distribution plots to be generated", nplots)
  if (nplots == 0) {
    df <- data.frame()
    dist.plots[[1]] <- ggplot(df) + geom_point() + xlim(0, 3) + ylim(0, 100) + ggtitle("No data to plot")
  }
  ## require(gridExtra)
  tmp <- do.call(marrangeGrob, c(dist.plots, list(nrow=nrow, ncol=ncol)))
  ggsave(filename = fn, tmp, width = 8, height = 16)
}


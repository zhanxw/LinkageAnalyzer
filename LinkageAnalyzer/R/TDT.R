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


#' Perform TDT
#' @keywords internal
TDT <- function(data, G2, gene, bin, log_file) {
  ############# continuous response doesn't work for TDT ###########

  if (bin == F) {
    return(1)
  }

  ############## clean data ##################

  # extract G2 genotype, mothers with NA value in genotype are discarded
  G2_gene <- unlist(G2[paste(gene$Gene, gene$Coordination), ])
  G2_gene <- G2_gene[-c(1:3)]
  G2_gene <- G2_gene[G2_gene %in% c("HET", "REF", "VAR")]

  # clean data
  data$mother <- as.vector(data$mother)
  data <- data[!is.na(data$gt), ]
  data <- data[!is.na(data$pt), ]

  ############# iTDT ####################

  # find gt of mothers
  data_iTDT <- data[data$mother %in% names(G2_gene), ]
  data_iTDT$mother_gt <- unlist(G2_gene[match(data_iTDT$mother, names(G2_gene))])
  data_iTDT$mother_gt <- convert_gt(data_iTDT$mother_gt, "additive")
  if (any(data_iTDT$mother_gt + 1 < data_iTDT$gt)) {
    report("w", paste("G2 genotype error for", gene[, 1]), log_file)
  }

  bc <- c()  # each element is the number of transmission informative for a (+) or
  # for A (-), from one parent

  for (mother in unique(data_iTDT$mother)) {
    ## print (mother)
    data_iTDT_mother <- data_iTDT[data_iTDT$mother == mother, ]
    if (sum(data_iTDT_mother$pt == "AFFECTED") == 0) {
      next
    }  # at least one affected

    q_A <- sum(data_iTDT_mother$pt == "AFFECTED" & data_iTDT_mother$gt == 0)
    q_U <- sum(data_iTDT_mother$pt != "AFFECTED" & data_iTDT_mother$gt == 0)
    p_A <- sum(data_iTDT_mother$pt == "AFFECTED" & data_iTDT_mother$gt == 1)
    p_U <- sum(data_iTDT_mother$pt != "AFFECTED" & data_iTDT_mother$gt == 1)
    r_A <- sum(data_iTDT_mother$pt == "AFFECTED" & data_iTDT_mother$gt == 2)
    r_U <- sum(data_iTDT_mother$pt != "AFFECTED" & data_iTDT_mother$gt == 2)

    if (data_iTDT_mother$mother_gt[1] == 0) {
      bc <- c(bc, q_A + p_U - p_A - q_U)  # F
    } else {
      bc <- c(bc, q_A + r_U - q_U - r_A)  # F
      bc <- c(bc, q_A + r_U - q_U - r_A)  # M
    }
  }

  bc <- bc[bc != 0]
  if (length(bc) > 0) {
    bc <- -bc
  }

  ############### sTDT #########################

  data_sTDT <- data[!data$mother %in% names(G2_gene), ]

  bcs <- c()  # each element is the difference in average affected allele count for each sibship

  for (mother in unique(data_sTDT$mother)) {
    data_sTDT_mother <- data_sTDT[data_sTDT$mother == mother, ]
    if (length(unique(data_sTDT_mother$pt)) == 1)
      {
        next
      }  # at least one affected and one unaffected

    mA <- sum(data_sTDT_mother$gt[data_sTDT_mother$pt == "AFFECTED"])/
        sum(data_sTDT_mother$pt == "AFFECTED")
    mU <- sum(data_sTDT_mother$gt[data_sTDT_mother$pt == "UNAFFECTED"])/
        sum(data_sTDT_mother$pt == "UNAFFECTED")
    bcs <- c(bcs, mA - mU)
  }

  bcs <- bcs[bcs != 0]
  if (length(bcs) > 0) {
    bcs <- bcs/abs(bcs)
  }

  ############### chisq-test ###################

  bc_combined <- c(bc, bcs)
  if (length(bc_combined) == 0) {
    return(1)
  }

  pval <- 1 - pchisq(sum(bc_combined)^2/sum(bc_combined^2), 1)
  return(convert_tail(sum(bc_combined) < 0, pval, "decreasing"))
}

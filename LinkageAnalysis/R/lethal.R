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


# G1 female is REF, G1 male is HET 1) Without G2 data G2 female (offspring) is
# REF/HET with 50% probability G2 male is HET (the same as G1 male)
# REF:HET:VAR=3:4:1 2) if G2 mother is REF, REF:HET:VAR=1:1:0 3) if G2 mother is
# HET, REF:HET:VAR=1:2:1


#' Calculate single lethal p-value
#' @param nvarFromHet counts of G3 genotyped and with HET mom
#' @param nFromHet     counts of G3 with VAR and with HET mom
#' @param nvarFromUnknown counts of offsprintgs genotyped and with Unknown mom
#' @param nFromUnknown     counts of G3 with VAR and with Unknown mom
single.lethal.getPvalue <- function(nvarFromHet, nFromHet,
                                    nvarFromUnknown, nFromUnknown) {
  #' g2Unknown can be ungenotyped or g2 = VAR
  #' @param nvarFromHet
  getP <- function(nvarFromHet, nFromHet, nvarFromUnknown, nFromUnknown) {
    ##dbinom(nvarFromHet, nFromHet, 0.25) * getP.g2unknown(nvarFromUnknown, nFromUnknown)
    dbinom(nvarFromHet, nFromHet, 0.25) * dbinom(nvarFromUnknown, nFromUnknown, 0.125)
  }

  nvar <- nvarFromHet + nvarFromUnknown
  p <- 0.
  for (i in seq(0, nvar)) {
    j <- seq(0, i)
    tmp <- getP(j, nFromHet, i - j , nFromUnknown)
    ## print(i)
    ## print(j)
    ## print(tmp)
    p <- p + sum(tmp)
  }
  p
}

# this function tests whether each gene is a homozygous lethal with G2 data the
# H0 is a gene is not homozygous lethal
single_lethal <- function(genotype, genes, phenotype, G2) {
  if (FALSE) {
    wd <- getwd()
    save(file = "single_lethal.Rdata", list = ls())
  }

  lethal <- rep(1, dim(genotype)[1])
  mothers <- as.vector(unique(phenotype$mother))

  for (i in 1:dim(genotype)[1]) {
    ## if (TRUE) {
    ##   cat("DBG:", "i = ", i, " ", genes$Gene[i], "\n")
    ## }

    # skip genes with too many invalid values or on chrX
    gt <- unlist(genotype[i, ])
    gt <- convert_gt(gt, "recessive")
    if (sum(!is.na(gt)) < 3 || genes$chr[i] == "X") {
      next
    }

    # aggregate counts from mother and its offsprings
    nVarFromHet <- 0
    nFromHet <- 0
    nVarFromUnknown <- 0
    nFromUnknown <- 0

    for (mother in mothers) {
      mother_gt <- G2[paste(genes$Gene[i], genes$Coordination[i]), mother]
      # in case G2 genotype is unknown, if any child is VAR, the mother is definitely
      # HET
      if (any(gt[phenotype$mother == mother & !is.na(gt)] == 2)) {
        cat("Fix mother genotype\n")
        mother_gt <- "HET"
      }

      if (is.null(mother_gt)) {
        nVarFromUnknown <- nVarFromUnknown + length(gt[phenotype$mother == mother & (gt==2)])
        nFromUnknown <- nFromUnknown + length(gt[phenotype$mother == mother & !is.na(gt)])
      } else if (mother_gt == "REF") {
        ## do nothing
      } else if (mother_gt == "HET") {
        nVarFromHet <- nVarFromHet + length(gt[phenotype$mother == mother & (gt==2)])
        nFromHet <- nFromHet + length(gt[phenotype$mother == mother & !is.na(gt)])
      } else if (mother_gt == "VAR") {
        ## should not happen
      } else {
        ## record unknown types
        cat("unexpected mother_gt = ", mother_gt, "\n")
        nVarFromUnknown <- nVarFromUnknown + length(gt[phenotype$mother == mother & (gt==2)])
        nFromUnknown <- nFromUnknown + length(gt[phenotype$mother == mother & !is.na(gt)])
      }
      ## cat("DBG", " mother = ", mother, " mother_gt = ", mother_gt, "\n")
      ## print(table(gt[phenotype$mother == mother], exclude =  NULL))
    }
    ## print(c(nVarFromHet, nFromHet, nVarFromUnknown, nFromUnknown))
    lethal[i] <- single.lethal.getPvalue(nVarFromHet, nFromHet, nVarFromUnknown, nFromUnknown)
  }
  return(lethal)
}

# this function runs random sampling to find the number of G3 mice with desired
# genotype
double_sample <- function(mother_gt1, mother_gt2, n_G3) {
  convertToGenotypeChar <- function(x) {
    if (is.na(x)) {
      return ("unknown")
    }
    if (is.character(x)) {
      return(x)
    }
    if (x > 1.5) { return("VAR")}
    if (x < 0.5) { return("REF")}
    return("HET")
  }
  mother_gt1 <- convertToGenotypeChar(mother_gt1)
  mother_gt2 <- convertToGenotypeChar(mother_gt2)

  if (is.null(mother_gt1) || is.na(mother_gt1) ||
      mother_gt1 %in% c("FALSE", "FAILED", "ERROR")) {
    mother_gt1 <- "unknown"
  }
  if (is.null(mother_gt2) || is.na(mother_gt2) ||
      mother_gt2 %in% c("FALSE", "FAILED", "ERROR")) {
    mother_gt2 <- "unknown"
  }

  prob <- -1
  if (mother_gt1 == "REF" && mother_gt2 == "REF") {
    # cannot give the genotype we desired for sure
    ##return(rep(0, n_trial))
    prob <- 0
  }
  if (mother_gt1 == "HET" && mother_gt2 == "HET") {
    prob <- 5/16
  }
  if ((mother_gt1 == "HET" && mother_gt2 == "REF") ||
      (mother_gt1 == "REF" && mother_gt2 == "HET")) {
    prob <- 1/8
  }
  if (mother_gt1 == "unknown" && mother_gt2 == "unknown") {
    prob <- (1/8 + 1/8 + 5/16)/4
  }
  if ((mother_gt1 == "unknown" && mother_gt2 == "REF") ||
      (mother_gt1 == "REF" && mother_gt2 == "unknown")) {
    prob <- (1/8)/2
  }
  if ((mother_gt1 == "unknown" && mother_gt2 == "HET") ||
      (mother_gt1 == "HET" && mother_gt2 == "unknown")) {
    prob <- (1/8 + 5/16)/2
  }

  ## return(rbinom(n = n_trial, size = n_G3, prob = prob))
  if (prob < 0) {
    snapshot("double_sample", "dbg.double_sample.Rdata", force = TRUE)
    logwarn("double_sample() encountered a strange case")
    return(0)
  }
  return(prob)
}

## calculate P(sum(X_i) <= obs)
## when X_i ~ binomial(size[i], p[i])
calculate.prob <- function (obs, size, p) {
  if (FALSE) {
    cat("DBG in calculate.prob()\n")
    wd <- getwd()
    save(list=ls(), file = "dbg.calculate.prob.Rsave")
  }
  N <- length (obs)
  if (N == 0) {
    return(1)
  }
  ret <- list ()
  n.obs <- sum (obs)
  for (i in 1:N) {
    ret [[i]] <- dbinom (0:n.obs, size = size [i], prob = p [i])
  }

  combine.prob <- function (v1, v2) {
    stopifnot (is.numeric (v1))
    stopifnot (is.numeric (v2))
    n <- length (v1)
    #print (v1)
    #print (v2)
    stopifnot (n == length (v2))
    ret <- rep (NA, n)
    for (i in 1:n) {
      ret [i] <- sum (v1 [1:i] * rev (v2 [1:i]))
    }
    ret
  }
  ## combine.prob (ret [[1]], ret [[2]])
  reduce.prob <- function (l) {
    stopifnot (is.list (l))
    n <- length (l)
    if (n < 2) {
      return (l [[1]])
    }
    v <- NULL
    for (i in 1: (n-1)) {
      if (i == 1) {
        v <- combine.prob (ret [[i]], ret [[i + 1]])
      } else {
        v <- combine.prob (v, ret [[i+1]])
      }
      #    cat (i, " : ")
      #    print (v)
    }
    return (v)
  }
  tmp <- reduce.prob (ret)
  sum (tmp)
}

# this function tests for synthetic lethality of two genes (VAR,HET; HET,VAR; and
# VAR,VAR)
double_lethal <- function(data, input, i, j) {
  mothers <- as.vector(unique(data$mother))
  if (any(input$genes[c(i, j), "chr"] == "X")) {
    return(1)
  }  # don't handle gene on chrX for the moment

  ## MonteCarlo <- matrix(data = 0, nrow = n_trial, ncol = length(mothers))
  ## colnames(MonteCarlo) <- mothers
  prob <- rep(0, length(mothers))
  numG3 <- rep(0, length(mothers))
  numObs <- rep(0, length(mothers))
  names(prob) <- names(numG3) <- names(numObs) <- mothers

  for (mother in mothers) {
    mother_gt1 <- input$G2[paste(input$genes$Gene[i], input$genes$Coordination[i]),
                           mother]
    mother_gt2 <- input$G2[paste(input$genes$Gene[j], input$genes$Coordination[j]),
                           mother]
    n_G3 <- sum(data$mother == mother)
    # observfed HET/VAR, VAR/HET or VAR/VAR
    tmp <- data[data$mother == mother, ]
    obs <- sum((tmp$gt1 + tmp$gt2) >= 3)

    # in case G2 genotype is unknown, if any child is VAR, the mother is definitely
    # HET
    if (any(data$gt1[data$mother == mother & !is.na(data$gt1)] == 2)) {
      mother_gt1 <- "HET"
    }
    if (any(data$gt2[data$mother == mother & !is.na(data$gt2)] == 2)) {
      mother_gt2 <- "HET"
    }
    ## MonteCarlo[, mother] <- double_sample(mother_gt1, mother_gt2, n_G3, n_trial)
    ## system.time(MonteCarlo[, mother] <- double_sample(mother_gt1, mother_gt2, n_G3, n_trial))
    prob[mother]  <- double_sample(mother_gt1, mother_gt2, n_G3)
    numG3[mother] <- n_G3
    numObs[mother] <- obs
  }

  ## collapsing
  tmp <- data.frame(prob = prob, numG3 = numG3, numObs = numObs)
  ## require(plyr)
  tmp <- ddply(tmp, .(prob), function(x) {c(numG3 = sum(x$numG3), numObs = sum(x$numObs))})
  tmp <- subset(tmp, prob != 0.0)
  pval <- calculate.prob(obs = tmp$numObs, size = tmp$numG3, p = tmp$prob)

  return(pval)
}

## calculate lehtal probability for double linkage
double.lethal.get.pvalue <- function(data) {
  mothers <- unique(data$mother)
  mothers <- mothers[mothers != "."]
  g2 <- data$iid[data$gen == 2]
  g2.mothers <- intersect(g2, mothers)
  if (length(g2.mothers) == 0) {
    logwarn("Cannot find any valid G2 mothers, skipped lethal calculation")
    return(0)
  }

  ## Sex chromsome not supported, and should be handled outside of this function
  ## if (any(input$genes[c(i, j), "chr"] == "X")) {
  ##   return(1)
  ## }  # don't handle gene on chrX for the moment

  ## MonteCarlo <- matrix(data = 0, nrow = n_trial, ncol = length(mothers))
  ## colnames(MonteCarlo) <- mothers
  prob <- rep(0, length(mothers))
  numG3 <- rep(0, length(mothers))
  numObs <- rep(0, length(mothers))
  names(prob) <- names(numG3) <- names(numObs) <- mothers

  for (mother in g2.mothers) {
    ## print(mother)
    mother_gt1 <- data$gt1[data$iid == mother]
    mother_gt2 <- data$gt2[data$iid == mother]
    if (length(mother_gt1) == 0) {
      logwarn("%s may not be a valid G2 mother", mother)
      next
    }
    n_G3 <- sum(data$mother == mother)
    # observfed HET/VAR, VAR/HET or VAR/VAR
    tmp <- data[data$mother == mother, ]
    obs <- sum((tmp$gt1 + tmp$gt2) >= 3)
    ## print(obs)
    if (is.na(obs)) {
      logwarn("No offsprings found from mother %s", mother)
      next
    }

    # the maternal-offspring genotypes should have been fixed, but add sanity check codes anyway
    # in case G2 genotype is unknown, if any child is VAR, the mother is definitely HET
    if (any(data$gt1[data$mother == mother & !is.na(data$gt1)] == 2)) {
      ## mother_gt1 <- "HET"
      if (is.na(mother_gt1) || (mother_gt1 < 0.5) || (mother_gt1 > 1.5)) {
        logwarn("mother genotype is not HET for %s", mother)
      }
    }
    if (any(data$gt2[data$mother == mother & !is.na(data$gt2)] == 2)) {
      ## mother_gt2 <- "HET"
      if (is.na(mother_gt2) || (mother_gt2 < 0.5) || (mother_gt2 > 1.5)) {
        logwarn("mother genotype is not HET for %s", mother)
      }
    }
    ## print(mother_gt1)
    ## print(mother_gt2)
    ## print(n_G3)
    ## print(prob[mother])
    ## print(double_sample(mother_gt1, mother_gt2, n_G3))
    prob[mother]  <- double_sample(mother_gt1, mother_gt2, n_G3)
    numG3[mother] <- n_G3
    numObs[mother] <- obs
  }

  ## collapsing
  tmp <- data.frame(prob = prob, numG3 = numG3, numObs = numObs)
  ## require(plyr)
  tmp <- ddply(tmp, .(prob), function(x) {c(numG3 = sum(x$numG3), numObs = sum(x$numObs))})
  tmp <- subset(tmp, prob != 0.0)
  pval <- calculate.prob(obs = tmp$numObs, size = tmp$numG3, p = tmp$prob)

  return(pval)
}

#' Examine mother-offspring pair and fix mother's genotype by her offsprings.
#' e.g. mom has VAR offsprings, mom genotype cannot be REF => mom must be HET
#' @param pheno ped data
#' @param geno.name specify genotype column headers
#' @keywords internal
crossCheckMotherGenotype <- function(pheno, geno.name = "gt") {
  mother <- NULL # bypass CRAN warning
  ## count mom's offspring phenotypes
  mom <- ddply(pheno, .(mother), function(x) {
    c(
        ref = sum(x[[geno.name]] == 0, na.rm = TRUE),
        var = sum(x[[geno.name]] > 1.5, na.rm = TRUE))
  })
  var <- NULL # bypass CRAN warning
  mom <- subset(mom, var > 0)
  total.fix <- 0
  fixed.mom.id <- vector("character", 0)
  for (i in seq_len(nrow(mom))) {
    mom.name <- mom[i, "mother"]

    if (!mom.name %in% pheno$iid ) {
      next
    }
    if (!is.na(pheno[ pheno$iid == mom.name, geno.name]) &&
        pheno[ pheno$iid == mom.name, geno.name] != 1) {
      pheno[ pheno$iid == mom.name, geno.name] <- 1 ## set to het
      loginfo("Fix mother %s by her offspring %s", mom.name, pheno$iid)
      total.fix <- total.fix + 1
      fixed.mom.id <- c(fixed.mom.id, mom.name)
    }
  }
  if (total.fix > 0) {
    loginfo("Fixed %d mother genotypes: %s", total.fix, paste0(unique(fixed.mom.id), collapse = ","))
  }
  pheno
}

countForLethal <- function(pheno) {
  nHetMom <- nVarFromHetMom <- nUnknownMom <- nVarFromUnknownMom <- 0
  mothers <- unique(pheno[, "mother"])
  g2 <- pheno$iid[pheno$gen == 2]
  g2.mothers <- intersect(g2, mothers)
  if (length(g2.mothers) == 0) {
    list(nVarFromHetMom, nHetMom, nVarFromUnknownMom, nUnknownMom)
  }
  ## loop each mother
  for (mother in g2.mothers) {
    if (mother == ".") { next } ## this should not happen
    idx <- which(mother == pheno$iid)
    if (length(idx) != 1) {
      ## mother not in ped
      ## nUnknownMom <- nUnknownMom + 1
      warning("Mother not in PED or multiple moms!")
      next
    }

    mother.gt <- pheno[idx, "gt"]
    if (is.na(mother.gt)) {
      offspring <- pheno[pheno$mother == mother & !is.na(pheno$gt), ]
      if (nrow(offspring) > 0 ) {
        nUnknownMom        <- nUnknownMom +        sum(!is.na(offspring$gt), na.rm = TRUE)
        nVarFromUnknownMom <- nVarFromUnknownMom + sum(offspring$gt == 2, na.rm = TRUE)
      }
      if (is.debug.mode()) {
        print(mother)
        print(offspring)
      }
    } else if (mother.gt == 0) {
      next
    } else if (mother.gt == 1) {
      offspring <- pheno[pheno$mother == mother & !is.na(pheno$gt), ]
      if (nrow(offspring) > 0 ) {
        nHetMom        <- nHetMom        + sum(!is.na(offspring$gt), na.rm = TRUE)
        nVarFromHetMom <- nVarFromHetMom + sum(offspring$gt == 2, na.rm = TRUE)
      }
      if (is.debug.mode()) {
        print(mother)
        print(offspring)
      }
    } else if (mother.gt == 2) {
      warning("Observe mother GT = HET :", mother)
      next ## should not happen ...
    } else {
      stop("Unrecognized mother genotype!")
    }
    ## print(mother)
    ## print(list(nHetMom, nVarFromHetMom, nUnknownMom, nVarFromUnknownMom))
  }
  list(nVarFromHetMom, nHetMom, nVarFromUnknownMom, nUnknownMom)
}

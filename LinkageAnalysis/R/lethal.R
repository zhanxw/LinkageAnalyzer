# G1 female is REF, G1 male is HET 1) Without G2 data G2 female (offspring) is
# REF/HET with 50% probability G2 male is HET (the same as G1 male)
# REF:HET:VAR=3:4:1 2) if G2 mother is REF, REF:HET:VAR=1:1:0 3) if G2 mother is
# HET, REF:HET:VAR=1:2:1


#' Calculate single lethal p-value
#' @param nvarFromHet counts of parents with HET
#' @param nvarHet     counts of offsprings with VAR and with HET mom
#' @param nvarFromUnknown counts of parents with Unknown genotype
#' @param nvarUnknown     counts of offsprings with VAR and with Unknown mom
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
  if (is.null(mother_gt1) || is.na(mother_gt1) ||
      mother_gt1 %in% c("FALSE", "FAILED", "ERROR")) {
    mother_gt1 <- "unknown"
  }
  if (is.null(mother_gt2) || is.na(mother_gt2) ||
      mother_gt2 %in% c("FALSE", "FAILED", "ERROR")) {
    mother_gt2 <- "unknown"
  }

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
  ret <- list ()
  n.obs <- sum (obs)
  for (i in 1:N) {
    ret [[i]] <- dbinom (0:n.obs, size = size [i], p = p [i])
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

  if (FALSE)  {
    cat ("DEBUG\n")
    wd <- getwd()
    fn <- "double_lethal.Rdata"
    save(list = ls(), file = fn)
    cat ( "double_lethal: data saved, ")
    cat (normalizePath(fn))
    cat ("\n")
  }
  if (FALSE) {
    load("/home/zhanxw/test.run/Rpackage.FACS_screen_B1b_cells.R0511.test1/double_lethal.Rdata")
  }
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
  library(plyr)
  tmp <- ddply(tmp, .(prob), function(x) {c(numG3 = sum(x$numG3), numObs = sum(x$numObs))})
  tmp <- subset(tmp, prob != 0.0)
  pval <- calculate.prob(obs = tmp$numObs, size = tmp$numG3, p = tmp$prob)

  ## obs <- sum((data$gt1 + data$gt2) >= 3)
  ## ## use a bit of adaptive
  ## mean <- sum(numG3 * prob)
  ## var <- sum(numG3 * prob * (1-prob))
  ## pval.approx <- pnorm( (obs - mean) /sqrt(var) )
  ## if (is.na(pval.approx)) {
  ##   #e.g. obs = mean = var = 0,
  ##   pval <- 1
  ## } else if (pval.approx > 0.1) {
  ##   pval <- doubleSample(prob, numG3, obs, 1000)
  ##   if (pval < 1e-4) {
  ##     ## this should rarely happen
  ##     pval <- doubleSample(prob, numG3, obs, 1000000)
  ##   }
  ## } else {
  ##   ## pval <- sum(apply(MonteCarlo, 1, sum) <= obs)/n_trial
  ##   pval <- doubleSample(prob, numG3, obs, 1000000)
  ## }
  return(pval)
}

##################################################
## OBSOLETE CODE
## # this function runs random sampling to find the number of G3 mice with VAR
## # genotype
## single_sample <- function(mother_gt, n_G3, n_trial) {
##   if (is.null(mother_gt) || is.na(mother_gt)) {
##     mother_gt <- "unknown"
##   }  # corresponding G2 data does not exist
##   if (mother_gt == "REF") {
##     return(0)
##   }
##   if (mother_gt == "HET") {
##     prob <- 1/4
##   }
##   if (!mother_gt %in% c("REF", "HET")) {
##     prob <- 1/8
##   }  # unknown, FAILED, FALSE
##   return(rbinom(n = n_trial, size = n_G3, prob = prob))
## }

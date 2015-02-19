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


#' Meta-analyze single gene in super pedigree
#'
#' Perform gene-based superpedigree analysis
#'
#' @param vcfFile genotype input file in VCF format
#' @param pedFile pedigree file in PED format
#' @param pheno.name phenotype name to analyze
#' @param output The output folder to put all the output files of the analysis
#' @param test The statistical test to be used for identifying significant genes
#' associated with phenotype. Default is "wG2", considering G2 mother
#' effect. Another option is "woG2", not considering that effect.
#' @param detect A character string specifying whether to detect clustering of
#' phenotypic scores and to transform continuous phenotypic scores into a binary
#' variable (affected and nonaffected). This parameter only works on continuous
#' phenotype scores. "never" (default): never transform; "always": always
#' transform; "auto": let the program decide whether to transform.
#' @param silent  Print intermediate messages to stdout if set to TRUE.
#' @param tail Either "decreasing", "increasing" or "both". "decreasing" tests
#' whether the mutation leads to decreased antibody reaction (default);
#' "increasing" tests whether the mutation leads to increased antibody reaction;
#' "both" tests the deviation in either direction.
#' @param prefix Default is "". An optional character string to be attached to
#' output file names. This can create personalized names for each job.
#' @param plot.it Default is TRUE, it controls whether to output linkage plots
#' and distribution plots.
#' @param transform.pheno Default is NULL. Use "log" if phenotypes need to log
#' transformed.
#' @param log.level Default WARN, but can be set from 'DEBUG', 'INFO', 'WARN'
#' and 'ERROR'.
#'
#' @return a list of three values will be returned:
#' the first one is returncode (0: success, non-zero: errors happened)
#' the second one is message, when return is non-zero, message gives explanations
#' the last one is analysis results.
#' @export
#' @examples
#' path <- system.file("extdata/gene",package="LinkageAnalyzer")
#' vcfFile <- file.path(path, "Snrnp40.mcmv.vcf")
#' pedFile <- file.path(path, "Snrnp40.mcmv.ped")
#' pheno.name <- "log"
#' output <- file.path(path, "output")
#' ret <- gene.single.link(vcfFile, pedFile, pheno.name, output,
#'                         test = "woG2", tail = "both", prefix="gene")
gene.single.link <- function(vcfFile, ## a vector of list
                             pedFile, ## a vector of list
                             pheno.name, ## which phenotype to use
                             output = ".",
                             test = "wG2",
                             detect = "never",
                             silent = T,
                             tail = "decreasing",
                             prefix = "",
                             plot.it = TRUE,
                             transform.pheno = NULL,
                             log.level = "WARN") {
  collectUsage("gene.single.link")
  log.file <- filename(output, prefix)$log_file
  ret <- tryCatch(
      {
        ret <- gene.single.link.impl(vcfFile, pedFile, pheno.name,
                                     output, test,
                                     detect, silent, tail, prefix,
                                     plot.it, transform.pheno,
                                     log.level)
        if (ret$returncode == 0) {
          msg <- paste("Exit successfully", ret$message, sep = " ")
        } else {
          if (isIgnorableError(ret)) {
            # this is a special error,
            # meaning we will treat it as normal exit but no output files
            # so returncode are changed from 1 to 0
            ret$returncode = 0
            msg <- paste("Exit successfully with no outputs due to", ret$message)
          } else {
            msg <- paste("Exit failed", ret$message, sep = " ")
          }
        }
        report("m", msg, log.file)
        return(ret)
      },
      error = function(err) {
        snapshot("gene.single.link.impl", "dbg.gene.single.link.impl.Rdata")

        if (err$message == "Response is constant - cannot fit the model"){
          # this is another special error,
          # meaning phenotypes are the same and statistical analysis cannot be performed
          # so returncode are changed from 1 to 0
          msg <- paste("Exit successfully with no outputs due to", err$message)
          report("m", msg, log.file)
          return(list(returncode = 0, message = msg))
        }

        reportError(err)
        msg <- ifelse(is.null(err[["message"]]), "UnknownError", err$message)
        msg <- paste("Exit failed", msg, sep = " ")
        report("m", msg, log.file)
        return(list(returncode = 1, message = msg))
      })

  return(ret)
}

gene.single.link.impl <- function(vcfFile, ## a vector of list
                                  pedFile, ## a vector of list
                                  pheno.name, ## which phenotype to use
                                  output = ".",
                                  test = "wG2",
                                  detect = "never",
                                  silent = T,
                                  tail = "decreasing",
                                  prefix = "",
                                  plot.it = TRUE,
                                  transform.pheno = NULL,
                                  log.level = 'WARN') {
  start.time <- Sys.time()

  ## set up log file
  log.file <- file.path(getwd(), file.path(output, "log.txt"))
  setupLogging(log.level, log.file)
  fns <- filename(output, prefix)  # generate output file names

  ## examine input parameters
  if (!tail %in% c("increasing", "decreasing", "both")) {
    msg <- paste0("ERROR: Unrecognized option for tail [ ", tail, "]\n")
    logerror(msg)
    return(list(returncode = 1, message = msg))
  }

  ## record running environment
  recordRunningInfo()
  
  ## read data
  tmp <- load.vcf.ped(vcfFile, pedFile, pheno.name)
  if (isSuccess(tmp)) {
    loginfo("VCF/PED loading succeed.")
    vcf <- tmp$vcf
    ped <- tmp$ped
  } else {
    loginfo("VCF/PED loading failed.")
    return(list(returncode = 1, message = msg))
  }
  snapshot("gene.link", "dbg.gene.Rdata")

  loginfo("Detecting phenotype type.\n")
  tmp <- process.phenotype(ped, pheno.name, detect)
  if (isSuccess(tmp)) {
    ped <- tmp$ped
    loginfo("Phenotype processed successfully.\n")
  } else {
    logerror("Phenotype has critical issues.\n")
    return(tmp)
  }

  ## make an output skeleton
  nSite <- length(vcf$CHROM)
  nas <- rep(NA, nSite)
  gene <- str_replace(sapply(vcf$INFO, function(x) {str_split(x, ";")[[1]][1]}), "GENE=", "")
  ret <- data.frame(Gene = gene, chr = vcf$CHROM, pos = vcf$POS,
                    REF = nas ,HET =nas, VAR = nas,
                    lethal = nas, additive = nas, recessive = nas, dominant = nas, TDT = nas,
                    Penetrance_REF = nas, Penetrance_HET = nas, Penetrance_VAR = nas, Semidominance = nas,
                    NMISS = nas, lethalCount = nas)
  rownames(ret) <- NULL


  ## calculate lethal
  type <- "lethal"
  loginfo("Perform %s test.", type)
  for (i in seq_len(nrow(ret))) {
    ## calculate lethal
    stopifnot(all(vcf$sampleId == ped$iid))
    if (i > 10 && is.debug.mode()) {
      loginfo("DEBUG skipped ", i, "th variant ..\n")
      next
    }
    loginfo("Access gene lethality for %s", gene[i])
    ## encode genotypes
    tmp <- ped
    tmp$gt <- convert_gt(vcf$GT[i,], "additive")
    ## check mother genotypes
    tmp <- crossCheckMotherGenotype(tmp)
    ## count mother, offspring by their genotypes
    tmp <- countForLethal(tmp)
    ret[i, "lethal"] <- do.call(single.lethal.getPvalue, tmp)
    ret[i, "lethalCount"] <- paste0(tmp, collapse = ":")
  }

  pheno.geno <- prepare.model.data(vcf, ped, pheno.name)
  pheno <- pheno.geno$pheno
  geno <- pheno.geno$geno
  nVariant <- nrow(geno)

  # set-up null model
  null <- create.null.model(pheno, pheno.name, test)
  if (!isSuccess(null)) {
    return(null)
  }
  has.random.effect <- null$has.random.effect
  isBinary <- null$isBinary

  snapshot("calc.genetic", "calc.genetic.Rdata")
  ret <- calc.genetic(ret, geno, pheno, pheno.name, isBinary)
  head(ret)

  # fit null model
  null <- fit.null.model(null)
  if (isSuccess(null)) {
    loginfo("Finished fitting null model\n")
  } else {
    logerror("Fitting null model failed!")
    return(null)
  }

  ## start fitting alternative models for each variant
  snapshot("gene.link", "dbg.gene.fit.alt.Rdata")

  dist.data <- list()
  dist.plots <- list()
  for (i in 1:nVariant) {
    if (i > 10 && is.debug.mode()) {
      cat("DEBUG skipped ", i, "th variant ..\n")
      next
    }
    loginfo("Process %d  out of %d variant: %s", i, nVariant, gene[i])

    # record data
    pheno$gt <- convert_gt(geno[i,], "additive")
    dist.data[[length(dist.data) + 1]] <- pheno

    ## skip mono site
    if (length(unique(pheno$gt)) == 1) {
      loginfo("Skip monomorphic site")
      next
    }

    # store null model
    null.orig <- null
    pheno.orig <- pheno

    # handle missing
    ret[i, "NMISS"] <- NMISS <- sum(is.na(pheno$gt))
    if (NMISS > 0 ) {
      loginfo("Refit null model due to %d missing genotypes", sum(is.na(pheno$gt)))
      pheno <- pheno[!is.na(pheno$gt), ]
      null <- fit.null.model(null.model, pheno)
      if (!isSuccess(null)) {
        logwarn("Null model fit failed.")
        null <- null.orig
        pheno <- pheno.orig
        next
      }
    }

    # calcualte each model
    raw.gt <- pheno$gt
    for (type in c("additive", "recessive", "dominant")) {
      loginfo("Perform %s test on %d samples", type, nrow(pheno))

      ## encode genotypes
      pheno$gt <- convert_gt(raw.gt, type)
      stopifnot(all(!is.na(pheno$gt)))

      ## skip this site when the genotypes are monomorphic
      if (length(unique(pheno$gt)[!is.na(unique(pheno$gt))]) <= 1) {
        loginfo("Skip monomorphic site under %s model", type)
        ret[i, type] <- pval <- 1
        next
      }

      ## fit alternative model
      alt.model <- paste0(null$formula, " + gt")
      if (has.random.effect) {
        alt <- run.mixed.effect.alt.model(isBinary, alt.model, pheno)
      } else {
        alt <- list(returncode = 1, message = "no random effect")
      }

      ## check convergence
      alt <- check.mixed.effect.alt.model.converge(alt)

      ## calculate p-value
      if (isSuccess(alt)) {
        ## no error occurred, using tradition anova tests
        pval <- anova(null$result, alt)$"Pr(>Chisq)"[2]
        if (is.na(pval)) {
          pval <- 1
        }
        direction <- fixef(alt)["gt"] > 0  # TRUE means protective effect, FALSE means harmful effect (desired)
        if (!is.na(direction)) {
          pval <- convert_tail(direction, pval, tail)
        }
      } else {
        loginfo("Refit alternative model using reduced data and Wald test\n")
        pval <- run.fixed.effect.alt.model(isBinary, alt.model, pheno, pheno.name, tail)
      }
      ret[i, type] <- pval
    }
    ret[i, "TDT"] <- NA     ## TODO: implement this

    ## restore null model
    null <- null.orig
    pheno <- pheno.orig

    # record graph
    if (any(ret[i, c("additive", "recessive", "dominant")] < 0.05, na.rm = TRUE)) {
      pheno$gt <- convert_gt(geno[i,], "additive")
      dist.title <- sprintf("%s (%s:%s)", gene[i], as.character(ret$chr[i]), as.character(ret$pos[i]))
      dist.plots[[length(dist.plots) + 1]] <- plot.distribution(pheno, pheno.name, dist.title)
    }
  }

  if (plot.it) {
    snapshot("plot.it", "plot.it.variant.Rdata")
    ## draw linkage plot
    linkage.plot.pdf <- changeSuffix(fns$linkage_file, ".pdf", ".variant.pdf")
    pdf(file = linkage.plot.pdf, height = 8, width = 20)
    plot.hist(pheno[,pheno.name], main = "Histogram of phenotype scores")
    plot.manhattan(data.frame(Chrom = ret$chr, Position = ret$pos, Gene = ret$Gene, Pval = ret$additive) , main = "Manhattan Plot: additive  model ")
    plot.manhattan(data.frame(Chrom = ret$chr, Position = ret$pos, Gene = ret$Gene, Pval = ret$recessive), main = "Manhattan Plot: recessive model ")
    plot.manhattan(data.frame(Chrom = ret$chr, Position = ret$pos, Gene = ret$Gene, Pval = ret$dominant) , main = "Manhattan Plot: dominant  model ")
    dev.off()
    loginfo(paste0("Generated ", linkage.plot.pdf))

    ## draw distribution plot
    dist.plot.pdf <- changeSuffix(fns$distrib_file, ".pdf", ".variant.pdf")
    save.dist.plot(dist.plots, nrow = 2, ncol = 1, fn = dist.plot.pdf)
    loginfo(paste0("Generated ", dist.plot.pdf))
  }

  ## Perform minimal p-value based approach
  loginfo("Summarize gene-level genetic mutation")
  snapshot("summarize.minP", "summarize.minP.Rdata")
  fam <- unique(pheno$fid)
  nFam <- length(fam)
  uniqGene <- unique(gene)
  genoCount <- matrix(NA, nrow = length(gene), ncol = nFam * 4)
  colName <- rep(fam, rep(4, nFam))
  colName <- paste(colName, rep(c("NA", 0, 1, 2), nFam), sep = "-")
  colnames(genoCount) <- colName
  for (i in seq_len(length(gene))) {
    for (j in seq_along(fam)) {
      f <- fam[j]
      idx <- pheno$fid == f
      g <- geno[i, idx]
      genoCount[i, j*4-3] <- sum(is.na(g))
      genoCount[i, j*4-2] <- sum(g == 0, na.rm = TRUE)
      genoCount[i, j*4-1] <- sum(g == 1, na.rm = TRUE)
      genoCount[i, j*4  ] <- sum(g == 2, na.rm = TRUE)
    }
  }
  genoCount <- cbind(data.frame(gene = gene), data.frame(genoCount))
  rownames(genoCount) <- NULL
  write.table(genoCount, file = changeSuffix(fns$csv_file, ".csv", ".preCollapseGeno.tbl"), quote = F, row.names = F)
  loginfo("Collapsing statistics")
  write.table(ret, file = changeSuffix(fns$csv_file, ".csv", ".variant.tbl"), quote = F, row.names = F)
  Gene <- NULL ## bypass CRAN check
  ret <- ddply(ret, .(Gene), function(x) {
    ret <- x[1,]
    ret$pos <-          natural.min(x$pos)
    ret$REF <-          natural.mean(x$REF)
    ret$HET <-          natural.mean(x$HET)
    ret$VAR <-          natural.mean(x$VAR)
    ret$lethal <-       natural.min(x$lethal)
    ret$additive <-     natural.min(x$additive)
    ret$recessive <-    natural.min(x$recessive)
    ret$dominant <-     natural.min(x$dominant)
    ret
  })
  rownames(ret) <- ret$Gene
  ret <- ret[uniqGene, ]
  tmp.fn <- changeSuffix(fns$csv_file, ".csv", ".minP.tbl")
  write.table(ret, file = tmp.fn, quote = F, row.names = F)
  loginfo(paste0("Generated ", tmp.fn))

  ## Loop each gene
  geno.collapse <- collapseGenotypeByGene(geno, gene)
  nGene <- nrow(geno.collapse)
  uniqGene <- rownames(geno.collapse)
  retGene <- ddply(ret, .(Gene), function(x) {
    data.frame(
        chr = unique(x$chr)[1],
        pos = min(as.integer(x$pos)),
        REF = NA,
        HET = NA,
        VAR = NA,
        lethal = natural.min(x$lethal), ## using minimal pval for lethal
        additive = NA,
        recessive = NA,
        dominant = NA,
        numVariant = nrow(x),
        variant = paste(x$pos, sep = ",")) })
  rownames(retGene) <- retGene$Gene
  retGene <- retGene[uniqGene, ]

  snapshot("calc.genetic.gene", "calc.genetic.gene.Rdata")
  retGene <- calc.genetic(retGene, geno.collapse, pheno, pheno.name, isBinary)
  head(retGene)

  ## start fitting alternative models for each variant
  snapshot("gene.link", "dbg.gene.fit.alt.gene.Rdata")

  dist.data <- list()
  dist.plots <- list()
  for (i in seq_len(nGene)) {
    if (i > 10 && is.debug.mode()) {
      cat("DEBUG skipped ", i, "th gene ..\n")
      next
    }
    loginfo("Process %d  out of %d gene: %s", i, nVariant, uniqGene[i])
    # record data
    pheno$gt <- geno.collapse[i,]
    dist.data[[length(dist.data) + 1]] <- pheno

    ## skip mono site
    if (length(unique(pheno$gt)) == 1) {
      loginfo("Skip monomorphic site")
      next
    }

    # store null model
    null.orig <- null
    pheno.orig <- pheno

    # handle missing
    retGene[i, "NMISS"] <- NMISS <- sum(is.na(pheno$gt))
    if (NMISS > 0 ) {
      loginfo("Refit null model due to %d missing genotypes", sum(is.na(pheno$gt)))
      pheno <- pheno[!is.na(pheno$gt), ]
      null <- fit.null.model(null.model, pheno)
      if (!isSuccess(null)) {
        logwarn("Null model fit failed.")
        null <- null.orig
        pheno <- pheno.orig
        next
      }
    }

    # calcualte each model
    raw.gt <- pheno$gt
    for (type in c("additive", "recessive", "dominant")) {
      loginfo("Perform %s test on %d samples", type, nrow(pheno))

      ## encode genotypes
      pheno$gt <- convert_gt(raw.gt, type)
      stopifnot(all(!is.na(pheno$gt)))

      ## skip this site when the genotypes are monomorphic
      if (length(unique(pheno$gt)[!is.na(unique(pheno$gt))]) <= 1) {
        loginfo("Skip monomorphic site under %s model", type)
        retGene[i, type] <- pval <- 1
        next
      }

      ## fit alternative model
      alt.model <- paste0(null$formula, " + gt")
      if (has.random.effect) {
        alt <- run.mixed.effect.alt.model(isBinary, alt.model, pheno)
      } else {
        alt <- list(returncode = 1, message = "no random effect")
      }

      ## check convergence
      alt <- check.mixed.effect.alt.model.converge(alt)

      ## calculate p-value
      if (isSuccess(alt)) {
        ## no error occurred, using tradition anova tests
        pval <- anova(null$result, alt)$"Pr(>Chisq)"[2]
        if (is.na(pval)) {
          pval <- 1
        }
        direction <- fixef(alt)["gt"] > 0  # TRUE means protective effect, FALSE means harmful effect (desired)
        if (!is.na(direction)) {
          pval <- convert_tail(direction, pval, tail)
        }
      } else {
        loginfo("Refit alternative model using reduced data and Wald test\n")
        pval <- run.fixed.effect.alt.model(isBinary, alt.model, pheno, pheno.name, tail)
      }
      retGene[i, type] <- pval
    }
    retGene[i, "TDT"] <- NA     ## TODO: implement this

    ## restore null model
    null <- null.orig
    pheno <- pheno.orig

    # record graph
    if (any(retGene[i, c("additive", "recessive", "dominant")] < 0.05, na.rm = TRUE)) {
      pheno$gt <- convert_gt(geno.collapse[i,], "additive")
      dist.title <- sprintf("%s (%s:%s)", gene[i], as.character(retGene$chr[i]), as.character(retGene$pos[i]))
      dist.plots[[length(dist.plots) + 1]] <- plot.distribution(pheno, pheno.name, dist.title)
    }
  }
  if (plot.it) {
    snapshot("plot.it.gene", "plot.it.gene.Rdata")
    ## draw linkage plot
    linkage.plot.pdf <- fns$linkage_file
    pdf(file = linkage.plot.pdf, height = 8, width = 20)
    plot.hist(pheno[,pheno.name], main = "Histogram of phenotype scores")
    plot.manhattan(data.frame(Chrom = retGene$chr, Position = retGene$pos, Gene = retGene$Gene, Pval = retGene$additive) , main = "Manhattan Plot: additive  model ")
    plot.manhattan(data.frame(Chrom = retGene$chr, Position = retGene$pos, Gene = retGene$Gene, Pval = retGene$recessive), main = "Manhattan Plot: recessive model ")
    plot.manhattan(data.frame(Chrom = retGene$chr, Position = retGene$pos, Gene = retGene$Gene, Pval = retGene$dominant) , main = "Manhattan Plot: dominant  model ")
    dev.off()
    loginfo(paste0("Generated ", linkage.plot.pdf))

    ## draw distribution plot
    dist.plot.pdf <- fns$distrib_file
    save.dist.plot(dist.plots, nrow = 2, ncol = 1, fn = dist.plot.pdf)
    loginfo(paste0("Generated ", dist.plot.pdf))
  }

  write.table(retGene, file = changeSuffix(fns$csv_file, ".csv", ".tbl"), quote = F, row.names = F)
  write.table(retGene, file = fns$csv_file, quote = F, row.names = F, sep = ",")
  loginfo(paste0("Generated ", fns$csv_file))


  wd <- getwd()
  save(list = ls(), file = fns$results_file)

  end.time <- Sys.time()
  diff.time <- difftime(end.time, start.time, units = "secs")
  msg <- (sprintf("gene.single.link.impl() finished in %.3f seconds - [vcfFile=%s;pedFile=%s;pheno.name=%s;test=%s;detect=%s;tail=%s]",
                  diff.time,
                  vcfFile,
                  pedFile,
                  pheno.name,
                  test,
                  detect,
                  tail))
  loginfo(msg)
  loginfo("Exit successfully")
  return(list(returncode = 0, message = "", result = ret))
}

## return collapsed genotype from @param geno, ordered by @param gene
collapseGenotypeByGene <- function(geno, gene) {
  snapshot("collapseGenotypeByGene", "collapseGenotypeByGene.Rdata")
  uniq.gene <- unique(gene)
  d <- data.frame(gene = gene, geno)
  ret <- daply(d, .(gene), function(x) {
    # print(x)
    apply(x[, -1, drop = FALSE], 2, function(x) {
      natural.max(x)
    })
  })
  ret <- ret[uniq.gene, ]
  ret
}

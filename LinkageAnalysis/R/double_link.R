#' @title Statistical analysis for testing the combinatory effect of two mutations
#'
#' @description This function tests whether the combination of two mutations
#' could be correlated with phenotype.
#' @param vcfFile genotype input file in VCF format
#' @param pedFile pedigree file in PED format
#' @param pheno.name phenotype name to analyze
#' @param output The output folder to put all the output files of the analysis
#' @param test The statistical test to be used for identifying significant genes
#' associated with phenotype. Default is "woG2", without considering G2 mother
#' effect. Another option is "wG2", considering that effect. "wG2" will take a
#' long time.
#' @param detect A character string specifying whether to detect clustering of
#' phenotypic scores and to transform continuous phenotypic scores into a binary
#' variable (affected and nonaffected). This parameter only works on continuous
#' phenotype scores. "never" (default): never transform; "always": always
#' transform; "auto": let the program decide whether to transform.
#' @param silent Print intermediate messages to stdout if set to TRUE.
#' @param tail Either "decreasing", "increasing" or "both". "decreasing" tests
#' whether the mutation leads to decreased antibody reaction (default);
#' "increasing" tests whether the mutation leads to increased antibody reaction;
#' "both" tests the deviation in either direction.
#' @param cutoff_single A numerical value smaller than 1. When between 0-1, only
#' genes whose single effect has a p value larger than \code{cutoff_single} will
#' be considered in testing for combinatory effect. Setting \code{cutoff_single}
#' to a negative number will turn off this screening (not recommended). Default:
#' 0.01.
#' @param prefix Default is "". An optional character string to be attached to
#' output file names. This can create personalized names for each job.
#' @param plot.it Default is TRUE, it controls whether to output linkage plots
#' and distribution plots.
#' @param transform.pheno Default is NULL. Use "log" if phenotypes need to log
#' transformed.
#' @param log.level Default WARN, but can be set from 'DEBUG', 'INFO', 'WARN'
#' and 'ERROR'.
#'
#' @details
#' The \code{double.link} function has 4 modes:
#' recessive, additive, dominant and inhibitory. Refer to the tutorial attached
#' in the package for the assumption in each mode. Recessive, additive, dominant
#' and inhibitory modes test whether the interaction between two genes could
#' explain the phenotypic variation with different assumptions of how the
#' interaction will take place.  Basically, ANOVA tests are applied with
#' different assumption in each mode to test the combinatory effects of two
#' genes. \code{cutoff_single} is suggested to be set to a value between 0-1, in
#' order to mask genes which themselves are already shown to be significantly
#' correlated with phenotype by \code{single.link}.  A PDF file
#' \code{linkage_plot} will be generated. It will contain 4 heatmap plots of p
#' values in each of the four modes. The p values are transformed to a -log10
#' scale. 4 CSV files containing the matrix of p values for each of the 4 modes
#' will be generated. A second PDF file \code{distribution_plot} will be
#' generated containing scattorplots of phenotypes scores or categorical
#' phenotypes for each pair of genes deemed to be significant. A TXT file
#' \code{log} containing messages, warnings and errors of the program will be
#' generated.  Synthetic lethality is tested in order to find whether G3 mice
#' with genotype of VAR,VAR; VAR,HET; and HET,VAR will have any decreased chance
#' of survival. The p values for synthetic lethality are presented along with
#' the p values for combinatory effect.
#' @return a list of two values will be returned. One is returncode (0: success).
#' The other is analysis results.
#' @seealso \code{\link{single.link}}.
#' @export
#' @examples
#' path <- system.file("extdata/double",package="LinkageAnalysis")
#' vcfFile <- file.path(path, "R0491_body_weight.vcf")
#' pedFile <- file.path(path, "R0491_body_weight.ped")
#' pheno.name <- "weight"
#' output <- file.path(path, "output")
#' ret <- double.link(vcfFile, pedFile, pheno.name, output,
#'                    test = "woG2", tail = "both", prefix="double")
double.link <- function(vcfFile,
                        pedFile,
                        pheno.name,
                        output = ".",
                        test = "woG2",
                        detect = "never",
                        silent = TRUE,
                        tail = "decreasing",
                        prefix = "",
                        cutoff_single = 0.01,
                        plot.it = TRUE,
                        transform.pheno = NULL,
                        log.level = 'WARN') {
  collectUsage("double.link")
  log.file <- filename(output, prefix)$log_file
  ret <- tryCatch(
      {
        ret <- double.link.impl(vcfFile, pedFile, pheno.name,
                                output, test, detect, silent, tail,
                                prefix, cutoff_single, plot.it, transform.pheno,
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
        snapshot("double.link.impl", "debug.double.link.impl.Rdata")
        reportError(err)
        msg <- ifelse(is.null(err[["message"]]), "UnknownError", err$message)
        msg <- paste("Exit failed", msg, sep = " ")
        report("m", msg, log.file)
        return(list(returncode = 1, message = msg))
      })

  return(ret)
}

double.link.impl <- function(vcfFile, pedFile, pheno.name,
                             output = ".", test = "woG2",
                             detect = "never",
                             silent = TRUE, tail = "decreasing", prefix = "",
                             cutoff_single = 0.01, plot.it = TRUE,
                             transform.pheno = NULL,
                             log.level = 'WARN') {
  start.time <- Sys.time()

  ## set up log file
  log.file <- file.path(getwd(), file.path(output, "log.txt"))
  basicConfig(log.level)
  addHandler(writeToFile, file = log.file)
  fns <- filename(output, prefix)  # generate output file names

  # check validity of parameters
  if (!is.numeric(cutoff_single) || cutoff_single > 1) {
    report("e", "Error in cutoff_single!", fns$log_file)
  }
  if (!tail %in% c("increasing", "decreasing", "both")) {
    report("e", "Unrecognized option for tail!", fns$log_file)
  }

  report("m", paste("Version:", packageVersion("LinkageAnalysis")), fns$log_file)
  report("m", paste("Date:", Sys.time()), fns$log_file)
  report("m", paste("Host:", Sys.info()["nodename"]) , fns$log_file)
  report("m", paste("Call:", deparse(sys.status()$sys.calls[[1]])), fns$log_file)


  ## read data
  snapshot("double.link.impl", "debug.double.before.load.Rdata")
  tmp <- load.vcf.ped(vcfFile, pedFile, pheno.name)
  if (isSuccess(tmp)) {
    loginfo("VCF/PED loading succeed.")
    vcf <- tmp$vcf
    ped <- tmp$ped
  } else {
    loginfo("VCF/PED loading failed.")
    return(list(returncode = 1, message = "VCF/PED loading failed."))
  }
  snapshot("double.link", "dbg.double.Rdata")

  nFamily <- length(unique(ped$fid))
  if (nFamily >  1) {
    msg <- sprintf("Detect %d families in PED file, please try other analysis: meta.single.link().", length(unique(ped$fid)))
    logerror(msg)
    return(list(returncode = 1, message = msg))
  } else if (nFamily <= 0){
    msg <- "No family detected!"
    logerror(msg)
    return(list(returncode = 1, message = msg))
  }

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
  signif <- matrix(data = 1, nrow = nSite, ncol = nSite)  # significance matrix
  colnames(signif) <- gene
  rownames(signif) <- gene
  ret <- list(recessive = signif, additive = signif,
              dominant = signif, inhibitory = signif,
              lethal = signif)

  ## ret <- data.frame(Gene = gene, chr = vcf$CHROM, pos = vcf$POS,
  ##                   REF = nas ,HET =nas, VAR = nas,
  ##                   lethal = nas, additive = nas, recessive = nas, dominant = nas, TDT = nas,
  ##                   Penetrance_REF = nas, Penetrance_HET = nas, Penetrance_VAR = nas, Semidominance = nas)
  ## rownames(ret) <- NULL

  loginfo("Load data complete")
  snapshot("double.link.impl", "debug.double.load.Rdata")


  ## calculate lethal (TODO)
  type <- "lethal"
  loginfo("Perform %s test.", type)
  for (i in 1:nSite) {
    if (is.debug.mode() && i > 10) {
      next
    }
    for (j in (i+1):nSite) {
      ## calculate lethal
      stopifnot(all(vcf$sampleId == ped$iid))
      ## skip sex chromosome
      if (any(!grepl(pattern = "(chr)?[0-9]+", vcf$CHROM[c(i,j)]))) {
        next
      }
      loginfo("Access gene lethality for %s x %s", gene[i], gene[j])
      ## encode genotypes
      tmp <- ped
      tmp$gt1 <- convert_gt(vcf$GT[i,], "additive")
      tmp$gt2 <- convert_gt(vcf$GT[j,], "additive")
      gen <- NULL ## bypass CRAN check
      tmp <- subset(tmp, gen == 2 | gen == 3)
      ## check mother genotypes
      tmp <- crossCheckMotherGenotype(tmp, "gt1")
      tmp <- crossCheckMotherGenotype(tmp, "gt2")
      ## ## count mother, offspring by their genotypes
      ## tmp <- countForDoubleLethal(tmp)
      # test for synthetic lethality
      double.lethal.get.pvalue(tmp)
      ret[["lethal"]][i, j] <- double.lethal.get.pvalue(tmp)
      ret[["lethal"]][j, i] <- ret[["lethal"]][i, j]
    }
  }

  pheno.geno <- prepare.model.data(vcf, ped, pheno.name)
  pheno <- pheno.geno$pheno
  geno <- pheno.geno$geno
  nVariant <- nrow(geno)

  # set-up null model
  null.model <- create.null.model(pheno, pheno.name, test)
  if (!isSuccess(null.model)) {
    return(null.model)
  }
  has.random.effect <- grepl("\\(", null.model)
  isBinary <- is.factor(pheno[,pheno.name])

  null <- fit.null.model(null.model, pheno, isBinary, has.random.effect)
  if (isSuccess(null)) {
    loginfo("Finished fitting null model\n")
  } else {
    logerror("Fitting null model failed!")
    loginfo("Outputting p value matrix anyway: %s",
            sub(pattern = "full", replacement = "\\*", fns$csv_file))
    for (type in c("recessive", "additive", "dominant", "inhibitory", "lethal")) {
      write.csv(ret[[type]], file = sub(pattern = "full", replacement = paste(type),
                                 fns$csv_file), quote = FALSE)
    }

    return(null)
  }

  ## start fitting alternative models for each variant
  snapshot("double.link", "dbg.double.fit.alt.Rdata")

  dist.data <- list()
  dist.plots <- list()
  # statistical test
  null.model <- NULL
  null.model.ok <- NULL
  report("m", paste("Total ", nSite, " gene(s) to test"), fns$log_file)
  for (i in 1:(nSite - 1)) {
    if (silent == FALSE) {
      if (is.debug.mode() && i > 10) {
        next
      }
      report("m", paste("---", gene[i], "---"), fns$log_file)
    }
    gt1 <- geno[i, ]  # genotype of the first gene

    for (j in (i + 1):nSite) {
      logdebug(sprintf("%s - %s x %s - (%d, %d, %d) - %.3f%%",
                    Sys.time(),
                    gene[i], gene[j],
                    i, j, nSite,
                    100. * ((j - i - 1) + (nSite - 1 + (nSite - i + 1)) * (i - 1) / 2) / (nSite * (nSite - 1)/2) ))
      gt2 <- geno[j, ]  # genotype of the second gene
      tmp <- convert_gt(gt1, "additive") - convert_gt(gt2, "additive")
      if (sum(!is.na(tmp)) <= 10) {
        logdebug("Small pedigree analysis starts.")
        ## next
      }  # too many NA values

      # prepare table of input and response variables
      data <- data.frame(pt = pheno[,pheno.name], sex = pheno$sex,
                         mother = pheno$mother, gt1 = convert_gt(gt1, "additive"),
                         gt2 = convert_gt(gt2, "additive"))
      hasMissingGeno <- any(is.na(data$gt1) | is.na(data$gt2))
      if (hasMissingGeno) {
        ## impute missing genotype to its mean
        naIdx <- is.na(data$gt1)
        data$gt1[naIdx] <- mean(data$gt1, na.rm = TRUE)
        naIdx <- is.na(data$gt2)
        data$gt2[naIdx] <- mean(data$gt2, na.rm = TRUE)
        ## data <- data[(!is.na(data$gt1)) & (!is.na(data$gt2)), ]
      }

      if (length(unique(data$gt1)) == 1 || length(unique(data$gt2)) == 1) {
        next
      }

      ## # if the two genes locate within 30MB, skipped
      ## if (input$genes$chr[i] == input$genes$chr[j] &&
      ##     abs(input$genes$pos[i] - input$genes$pos[j]) < 3e+07) {
      ##   next
      ## }

      # fit null model
      if (is.null(null.model)) {
        cat("fit null model\n")
        assign("last.warning", NULL, envir = baseenv())
        null.model <- anova_test(data, isBinary, test, silent, fns$log_file, tail, fit.null = TRUE)$null.model
        if (exists("last.warning", envir = baseenv()) &&
            !is.null(get("last.warning", baseenv()))) {
          report("w", "Null cannot be fitted!!", fns$log_file)
          null.model.ok <- FALSE
        } else {
          null.model.ok <- TRUE
        }
        cat("fit null model complete\n")
      }

      # test for combinatory effect
      for (type in c("recessive", "additive", "dominant", "inhibitory")) {
        if (type == "recessive") {
          data$gt <- (data$gt1 + data$gt2 >= 4) * 1
        }  else if (type == "additive") {
          data$gt <- data$gt1 * data$gt2
        } else if (type == "dominant") {
          data$gt <- (data$gt1 * data$gt2 >= 1) * 1
        } else if (type == "inhibitory") {
          data$gt <- data$gt1 * (data$gt2 == 0) + data$gt2 * (data$gt1 == 0)
        }
        if (length(unique(data$gt)) == 1) {
          next
        }  # no difference in predictor values
        if (FALSE) {
          system.time(anova_test(data, isBinary, test, silent, fns$log_file, tail, null.model = null.model))
          system.time(anova_test(data, isBinary, test, silent, fns$log_file, tail))
        }
        if (!null.model.ok)  {
          pval  <- 1
        } else {
          ## use tryCatch to avoid crashing
          pval <- tryCatch(
              {
                pval <- anova_test(data, isBinary, test, silent, fns$log_file, tail, null.model = null.model)$pvalue
              },
              error = function(err) {
                snapshot("double.link.impl", "debug.double.link.impl.Rdata")
                reportError(err)
                msg <- ifelse(is.null(err[["message"]]), "UnknownError", err$message)
                msg <- paste("Fitting failed", msg, sep = " ")
                report("m", msg, fns$log_file)
                return(list(returncode = 1, message = msg))
              })
          if (is.list(pval) && pval$returncode == 1) {
            ## error occured
            pval <- 1
          }

          ## cat('pval = ', pval , "\n")
          ret[[type]][i, j] <- pval
          ret[[type]][j, i] <- ret[[type]][i, j]
        }
      }

    } ## end loop j
  } ## end loop i

  ## # run single.link to determine which genes to mask
  ## if (silent == FALSE) {
  ##   report("m", "Running double.link", fns$log_file)
  ## }
  ## single.link(main_file = main_file, G2_file = G2_file, detect = detect, output = tempdir(),
  ##             silent = TRUE, test = test, tail = tail, plot.it = FALSE)  # run single.link
  ## SignifOfSingle <- read.csv(filename(tempdir(), prefix = "")$csv_file)  # get significance
  ## SignifOfSingle$inhibitory <- 1  # dummy for inhibitory mode, include all
  ## unlink(as.vector(filename(tempdir(), prefix = "")))
  SignifOfSingle <- list(recessive = rep(1, nSite),
                         additive = rep(1, nSite),
                         dominant = rep(1, nSite),
                         inhibitory = rep(1, nSite),
                         lethal = rep(1, nSite))

  # calculate bonferroni correction cutoff
  cutoff_pair <- 0.05/(nSite * (nSite - 1)/2)  # the cutoff is approximate using this rough estimation
  report("m", paste("The Bonferroni cutoff is", pretty_num(cutoff_pair)), fns$log_file)

  if (plot.it) {
    snapshot("double.plot.it", "dbg.double.plot.Rdata")
    # display results in heatmap (only genes that are not masked)
    if (silent == FALSE) {
      report("m", "Drawing heatmap", fns$log_file)
    }
    pdf(file = fns$linkage_file)
    par(mfrow = c(1, 1)) ## par() is necessary to draw multi-page heatmaps
    for (type in c("recessive", "additive", "dominant", "inhibitory", "lethal")) {
      # draw a heatmap of results
      ## mask <- SignifOfSingle[, type] > cutoff_single
      ## (TODO) need closer comparison with Tao's format
      ## heatmap_pval(sig, type, test, cutoff_pair, mask, fns)
      ## tmp <- data.frame(-log10(sig[[type]]))
      ## colnames(tmp) <- rownames(tmp) <- gene
      ## ggplot(tmp, aes(x = ))
      data <- as.matrix(-log10(ret[[type]]))
      calc.breaks <- function(data) {
        nSite <- nrow(data)
        if (length(data) == 0) {
          return(c(0, 1))
        }
        cutoff_pair <- 0.05/(nSite * (nSite - 1)/2)
        # break points for color scheme of heatmap
        rg <- range(data)
        if (rg[2] == Inf) {
          rg[2] <- -log10(cutoff_pair)
        }
        mid <- round(min(-log10(cutoff_pair), rg[2]/2))

        # the significant pairs should have better resolution
        a <- seq(from = 0, to = mid, by = mid/2)
        b <- seq(from = mid, to = max(ceiling(rg[2]), -log10(cutoff_pair)), length.out = 10)
        c <- unique(c(a, b))
        c <- c[order(c)]
        c
      }
      if (all(data == 0) || length(unique(as.vector(data))) == 1) {
        loginfo("Skip draw %d by %d heatmap for %s model as the p-values are constant.", dim(data)[1], dim(data)[2], type)
      } else {
        loginfo("Draw %d by %d heatmap for %s model", dim(data)[1], dim(data)[2], type)
        heatmap.2(data, dendrogram = "none", trace = "none",
                  Rowv = F, Colv = F, cexRow = 0.7, cexCol = 0.7,
                  main = paste("Heatmap of P values\n (", type, test, ")", sep = " "),
                  lmat = rbind(c(2, 3, 4), c(0, 1, 1)),
                  lhei = c(1, 5), lwid = c(0.5, 3, 1), keysize = 0.5,
                  density.info = "none",
                  breaks = calc.breaks(data))
      }
    }
    dev.off()
    loginfo(paste0("Generated ", fns$linkage_file))

    # make distribution plots (only genes that are not masked)
    if (silent == FALSE) {
      report("m", "Drawing distribution plots", fns$log_file)
    }
    pdf(file = fns$distrib_file)
    par(mfrow = c(3, 2))

    for (i in 1:(nSite - 1)) {
      for (j in (i + 1):nSite) {
        # whether to draw the distribution plot, which types pass the cutoff
        for (type in c("recessive", "additive", "dominant", "inhibitory", "lethal")) {
          ## Skip these masks
          ## if (ret[[type]][i, j] < cutoff_pair &&
          ##     SignifOfSingle[i, type] > cutoff_single &&
          ##     SignifOfSingle[j, type] > cutoff_single) {
          ##   flag <- paste(flag, " ", substr(type, 1, 2), " (", pretty_num(ret[[type]][i, j]), ")", sep = "")
          ## }
          if (ret[[type]][i, j] < 0.05) {
            gt1 <- geno[i, ]  # genotype of the first gene
            gt2 <- geno[j, ]  # genotype of the second gene
            data <- data.frame(pt = pheno[,pheno.name], sex = pheno$sex,
                               mother = pheno$mother, gt1 = convert_gt(gt1, "additive"),
                               gt2 = convert_gt(gt2, "additive"))
            data <- data[(!is.na(data$gt1)) & (!is.na(data$gt2)), ]

            # distribution of genotype vs. phenotype
            plot.double.distribution(data, isBinary, type, gene[c(i, j)])
          }
        }
      }
      ## plot(1:5, 1:5, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
    }
    dev.off()
  }

  # write p val matrix to file (full matrix)
  if (silent == FALSE) {
    report("m", "Outputting p value matrix", fns$log_file)
  }
  for (type in c("recessive", "additive", "dominant", "inhibitory", "lethal")) {
    write.csv(ret[[type]], file = sub(pattern = "full", replacement = paste(type),
                               fns$csv_file), quote = FALSE)
  }
  end.time <- Sys.time()
  diff.time <- difftime(end.time, start.time, units = "secs")
  msg <- (sprintf("double.link() finished in %.3f seconds - [vcfFile=%s;pedFile=%s;pheno.name=%s;test=%s;detect=%s;tail=%s]",
                  diff.time,
                  vcfFile, pedFile, pheno.name,
                  test,
                  detect,
                  tail))
  print(msg)
  report("m", msg, fns$log_file)

  return(list(returncode = 0, result = ret))
}

## double.link.impl.csv <- function(main_file, G2_file = "", output = ".", test = "woG2",
##                                  detect = "never",
##                                  silent = TRUE, tail = "decreasing", prefix = "",
##                                  cutoff_single = 0.01, plot.it = TRUE,
##                                  transform.pheno = NULL) {
##   start.time <- Sys.time()
##   # input
##   fns <- filename(output, prefix)  # generate output file names

##   report("m", paste("Version:", packageVersion("LinkageAnalysis")), fns$log_file)
##   report("m", paste("Date:", Sys.time()), fns$log_file)
##   report("m", paste("Host:", Sys.info()["nodename"]) , fns$log_file)
##   report("m", paste("Call:", deparse(sys.status()$sys.calls[[1]])), fns$log_file)

##   # read data,G2 is null is G2 dam genotype data are not available
##   snapshot("double.link.impl", "dbg.before.load.Rdata")
##   tmp <- get_data(main_file, G2_file, fns$log_file, detect, transform.pheno)
##   ## Before data load can be possibily fail, just write pval = 1 as results
##   if (!is.null(tmp$data$gene) ) {
##     ones <- rep(1, nrow(tmp$data$genes))
##     nas <- rep(NA, nrow(tmp$data$genes))

##     result <- tmp$data$genes
##     result$chr <- sub(pattern = "_.*", replacement = "", tmp$data$genes$Coordination, perl = TRUE)  # split Coordinates
##     result$pos <- sub(pattern = ".*_", replacement = "", tmp$data$genes$Coordination, perl = TRUE)

##     result$REF <- result$HET <- result$VAR <- nas
##     result$lethal <- result$additive <- result$recessive <- result$dominant <- result$TDT <- ones
##     result$Penetrance_REF <- result$Penetrance_HET <- result$Penetrance_VAR <- result$Semidominance <- nas
##     write.table(result, file = fns$csv_file, quote = F, row.names = F, sep = ",")
##   }
##   ## If data fails to load, quit
##   if (!is.null(ncol(tmp$data$genes) - 3)) {
##     nGeno <- ncol(tmp$data$genes) - 3
##   } else {
##     nGeno <- 0
##   }
##   if (tmp$returncode || nGeno <= 1) {
##     ## Data load failed, just write pval = 1 as results and quit
##     tmp$returncode = 0
##     return(tmp)
##   }

##   input <- tmp$data
##   report("m", "Load data complete", fns$log_file)

##   # check validity of parameters
##   if (!is.numeric(cutoff_single) || cutoff_single > 1) {
##     report("e", "Error in cutoff_single!", fns$log_file)
##   }
##   if (!tail %in% c("increasing", "decreasing", "both")) {
##     report("e", "Unrecognized option for tail!", fns$log_file)
##   }

##   # initialize significance matrix
##   snapshot("double.link.impl", "dbg.before.sig.Rdata")
##   signif <- matrix(data = 1, nrow = input$n, ncol = input$n)  # significance matrix
##   colnames(signif) <- input$genes$Gene
##   rownames(signif) <- input$genes$Gene
##   sig <- list(recessive = signif, additive = signif, dominant = signif,
##               inhibitory = signif,
##               lethal = signif)

##   # statistical test
##   null.model <- NULL
##   null.model.ok <- NULL
##   report("m", paste("Total ", input$n, " gene(s) to test"), fns$log_file)
##   for (i in 1:(input$n - 1)) {
##     ## if (i != 1) {
##     ##   cat("skip i = ", i , " in debug\n")
##     ##   next
##     ## }
##     if (silent == FALSE) {
##       report("m", paste("---", input$genes$Gene[i], "---"), fns$log_file)
##     }
##     gt1 <- unlist(input$genotype[i, ])  # genotype of the first gene

##     for (j in (i + 1):input$n) {
##       ## if (i > 10 || j > 10) {
##       ##   cat ("DBG skip\n")
##       ##   next
##       ## }
##       print(sprintf("%s - %s x %s - (%d, %d, %d) - %.3f%%",
##                     Sys.time(),
##                     input$genes$Gene[i], input$genes$Gene[j],
##                     i, j, input$n,
##                     100. * ((j - i - 1) + (input$n - 1 + (input$n - i + 1)) * (i - 1) / 2) /
##                     (input$n * (input$n - 1)/2) ))
##       gt2 <- unlist(input$genotype[j, ])  # genotype of the second gene
##       tmp <- convert_gt(gt1, "additive") - convert_gt(gt2, "additive")
##       if (sum(!is.na(tmp)) <= 10) {
##         next
##       }  # too many NA values

##       # prepare table of input and response variables
##       data <- data.frame(pt = input$phenotype$phenotype, sex = input$phenotype$sex,
##                          mother = input$phenotype$mother, gt1 = convert_gt(gt1, "additive"),
##                          gt2 = convert_gt(gt2, "additive"))
##       hasMissingGeno <- any(is.na(data$gt1) | is.na(data$gt2))
##       if (hasMissingGeno) {
##         ## impute missing genotype to its mean
##         naIdx <- is.na(data$gt1)
##         data$gt1[naIdx] <- mean(data$gt1, na.rm = TRUE)
##         naIdx <- is.na(data$gt2)
##         data$gt2[naIdx] <- mean(data$gt2, na.rm = TRUE)
##         ## data <- data[(!is.na(data$gt1)) & (!is.na(data$gt2)), ]
##       }

##       if (length(unique(data$gt1)) == 1 || length(unique(data$gt2)) == 1) {
##         next
##       }

##       # if the two genes locate within 30MB, skipped
##       if (input$genes$chr[i] == input$genes$chr[j] &&
##           abs(input$genes$pos[i] - input$genes$pos[j]) < 3e+07) {
##         next
##       }

##       # fit null model
##       ## source("/home/zhanxw/test.run/LinkageAnalysis/R/anova_test.R")
##       if (is.null(null.model)) {
##         cat("fit null model\n")
##         assign("last.warning", NULL, envir = baseenv())
##         null.model <- anova_test(data, input$bin, test, silent, fns$log_file, tail, fit.null = TRUE)$null.model
##         if (exists("last.warning", envir = baseenv()) &&
##             !is.null(get("last.warning", baseenv()))) {
##           report("w", "Null cannot be fitted!!", fns$log_file)
##           null.model.ok <- FALSE
##         } else {
##           null.model.ok <- TRUE
##         }
##         cat("fit null model complete\n")
##       }

##       # test for combinatory effect
##       for (type in c("recessive", "additive", "dominant", "inhibitory")) {
##         if (type == "recessive") {
##           data$gt <- (data$gt1 + data$gt2 >= 4) * 1
##         }  else if (type == "additive") {
##           data$gt <- data$gt1 * data$gt2
##         } else if (type == "dominant") {
##           data$gt <- (data$gt1 * data$gt2 >= 1) * 1
##         } else if (type == "inhibitory") {
##           data$gt <- data$gt1 * (data$gt2 == 0) + data$gt2 * (data$gt1 == 0)
##         }
##         if (length(unique(data$gt)) == 1) {
##           next
##         }  # no difference in predictor values
##         if (FALSE) {
##           system.time(anova_test(data, input$bin, test, silent, fns$log_file, tail, null.model = null.model))
##           system.time(anova_test(data, input$bin, test, silent, fns$log_file, tail))
##         }
##         if (!null.model.ok)  {
##           pval  <- 1
##           # cat("null.ok = ", null.model.ok, "\n")
##           ## cat('pval = ', pval , "\n")
##         } else {
##           ## use tryCatch to avoid crashing
##           pval <- tryCatch(
##               {
##                 pval <- anova_test(data, input$bin, test, silent, fns$log_file, tail, null.model = null.model)$pvalue
##               },
##               error = function(err) {
##                 snapshot("single.link.impl", "debug.single.link.impl.Rdata")
##                 print(str(err))
##                 print(err)
##                 msg <- ifelse(is.null(err[["message"]]), "UnknownError", err$message)
##                 msg <- paste("Fitting failed", msg, sep = " ")
##                 report("m", msg, fns$log_file)
##                 return(list(returncode = 1, message = msg))
##               })
##           if (is.list(pval) && pval$returncode == 1) {
##             ## error occured
##             pval <- 1
##           }

##           ## cat('pval = ', pval , "\n")
##           sig[[type]][i, j] <- pval
##           sig[[type]][j, i] <- sig[[type]][i, j]
##         }
##       }

##       # test for synthetic lethality
##       sig[["lethal"]][i, j] <- double_lethal(data, input, i, j)
##       sig[["lethal"]][j, i] <- sig[["lethal"]][i, j]
##     } ## end loop j
##   } ## end loop i

##   # run single.link to determine which genes to mask
##   if (silent == FALSE) {
##     report("m", "Running double.link", fns$log_file)
##   }
##   single.link(main_file = main_file, G2_file = G2_file, detect = detect, output = tempdir(),
##               silent = TRUE, test = test, tail = tail, plot.it = FALSE)  # run single.link
##   SignifOfSingle <- read.csv(filename(tempdir(), prefix = "")$csv_file)  # get significance
##   SignifOfSingle$inhibitory <- 1  # dummy for inhibitory mode, include all
##   unlink(as.vector(filename(tempdir(), prefix = "")))

##   # calculate bonferroni correction cutoff
##   cutoff_pair <- 0.05/(input$n * (input$n - 1)/2)  # the cutoff is approximate using this rough estimation
##   report("m", paste("The Bonferroni cutoff is", pretty_num(cutoff_pair)), fns$log_file)

##   if (plot.it) {
##     # display results in heatmap (only genes that are not masked)
##     if (silent == FALSE) {
##       report("m", "Drawing heatmap", fns$log_file)
##     }
##     ## par(mfrow = c(1, 1))
##     pdf(file = fns$pdf_file)

##     for (type in c("recessive", "additive", "dominant", "inhibitory", "lethal")) {
##       # draw a heatmap of results
##       mask <- SignifOfSingle[, type] > cutoff_single
##       heatmap_pval(sig, type, test, cutoff_pair, mask, fns)
##     }

##     dev.off()

##     # make distribution plots (only genes that are not masked)
##     if (silent == FALSE) {
##       report("m", "Drawing distribution plots", fns$log_file)
##     }
##     ## if (TRUE) {
##     ##   wd <- getwd()
##     ##   cat("DBG in ", wd, "\n")
##     ##   save(list = ls(), file = "tmp.Rdata")

##     ##   load("tmp.Rdata" ,verbose = TRUE)
##     ##   cutoff_pair <- 0.9
##     ##   sig[["recessive"]][1,2] <- 1e-14
##     ##   sig[["dominant"]][1,2] <- 1e-14
##     ##   sig[["inhibitory"]][1,2] <- 1e-14
##     ## }
##     pdf(file = fns$distrib_file)
##     par(mfrow = c(3, 2))

##     for (i in 1:(input$n - 1)) {
##       for (j in (i + 1):input$n) {
##         # whether to draw the distribution plot, which types pass the cutoff
##         flag <- ""
##         type.min <- NULL
##         p.min <- NULL
##         for (type in c("recessive", "additive", "dominant", "inhibitory", "lethal")) {
##           if (sig[[type]][i, j] < cutoff_pair &&
##               SignifOfSingle[i, type] > cutoff_single &&
##               SignifOfSingle[j, type] > cutoff_single) {
##             flag <- paste(flag, " ", substr(type, 1, 2), " (", pretty_num(sig[[type]][i, j]), ")", sep = "")
##           }
##           if (is.null(type.min) || p.min > sig[[type]][i, j]) {
##             type.min <- type
##             p.min <- sig[[type]][i, j]
##           }
##         }
##         if (flag == "") {
##           next
##         }
##         ## cat ("DBG flag = ", flag, "type.min", type.min, "\n")

##         gt1 <- unlist(input$genotype[i, ])  # genotype of the first gene
##         gt2 <- unlist(input$genotype[j, ])  # genotype of the second gene
##         data <- data.frame(pt = input$phenotype$phenotype, sex = input$phenotype$sex,
##                            mother = input$phenotype$mother, gt1 = convert_gt(gt1, "additive"),
##                            gt2 = convert_gt(gt2, "additive"))
##         data <- data[(!is.na(data$gt1)) & (!is.na(data$gt2)), ]

##         # distribution of genotype vs. phenotype
##         ## cat("DBG: draw double_distrib for type ", type, "\n")

##         double_distrib(flag, data, input, type.min, input$genes[c(i, j), ])
##       }
##     }

##     ## if (FALSE) {
##     ##   pdf(file = fns$distrib_file)
##     ##   par(mfrow = c(3, 2))
##     ##   double_distrib(flag, data, input, type, input$genes[c(i, j), ])  # distribution of genotype vs. phenotype
##     ##   plot(1:5, 1:5, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
##     ##   dev.off()
##     ## }
##     plot(1:5, 1:5, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
##     dev.off()
##     ## par(mfrow = c(1, 1))
##   }

##   # write p val matrix to file (full matrix)
##   if (silent == FALSE) {
##     report("m", "Outputting p value matrix", fns$log_file)
##   }
##   for (type in c("recessive", "additive", "dominant", "inhibitory", "lethal")) {
##     write.csv(sig[[type]], file = sub(pattern = "full", replacement = paste(type),
##                                fns$csv_file), quote = FALSE)
##   }
##   end.time <- Sys.time()
##   diff.time <- difftime(end.time, start.time, units = "secs")
##   msg <- (sprintf("double.link() finished in %.3f seconds - [main=%s;G2=%s;test=%s;detect=%s;tail=%s]",
##                   diff.time,
##                   main_file,
##                   G2_file,
##                   test,
##                   detect,
##                   tail))
##   print(msg)
##   report("m", msg, fns$log_file)

##   return(list(returncode = 0, result = sig))
## }

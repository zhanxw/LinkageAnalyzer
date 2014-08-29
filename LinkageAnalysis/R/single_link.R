#' @title Statistical analysis for correlation with phenotype and lethality on a single locus
#'
#' @description
#' This function reads the raw data and conducts statistical
#' analysis for the significance of correlation with phenotype and homozygous
#' lethality for each gene.
#'
#' @param main_file \code{raw_file} is an csv file containing main data.
#' See \code{details} for more information.
#' @param G2_file \code{G2_file} is an optional csv file containing G2 genotype
#' data, if such data are available. See \code{details} for more information.
#' @param output The output folder to put all the output files of the analysis
#' @param test The statistical test to be used for identifying significant genes
#' associated with phenotype. Default is "wG2", considering G2 mother effect.
#' Another option is "woG2", not considering that effect.
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
#'
#' @details First of all, this package depends on \code{lme4}. The current
#' version of this package on CRAN does not seem to be a stable one (at least on
#' my computer). I put an archived version of the lme4 package in the
#' \code{inst/extdata} folder of the \code{LinkageAnalysis} package. It works
#' well on my Linux system. You can install it from source on your Linux
#' platform, too. If you are working on PC/Mac, then you may also need to search
#' for one version of the \code{lme4} package that is stable on your machine.
#' \code{main_file} is a csv file with fixed format. The six cells in A1-B3 can
#' be anything. A4 must be "Gene", B4 must be "Coordination", C4 must be
#' "Amplicon", C3 must be "Phenotype", C2 must be "G3 gender" and C1 must be "G2
#' Dam eartag". The cells below A4 are gene names. The cells below B4 are
#' genomic locations ("1_100", for example, means chr1 and the 100th bp). The
#' cells below C4 are Amplicon. The cells to the right of C4 G3 eartag. The
#' cells to the right of C3 are phenotypes. They can be categorical values like,
#' "Affected", "Intermediate" and "Unaffected" or they can be numerical values,
#' but can not be a mixture of both. The cells to the right of C2 are G3 gender
#' and must be one of "F" and "M". The cells to the right of C1 are "G2 Dam
#' eartag". The cells in the bottom-right area are the genotypes. "REF", "HET"
#' and "VAR" correspond to wild-type, heterozygous mutant and homozygous
#' mutant. Other symbols will be treated as NA values. A sample raw input file
#' is stored in the \code{inst/extdata} folder of this package.  \code{G2_file}
#' is an optional csv file with fixed format. A1 must be "Gene", B1 must be
#' "Coordination" and C1 must be "Amplicon". The cells below A1 are gene
#' names. The cells below B1 are genomic locations. The cells below C1 are
#' Amplicon. The cells to the right of C1 are G2 dam names. The cells below G2
#' dam names are the G2 genotypes. "REF", "HET" and "VAR" correspond to
#' wild-type, heterozygous mutant and homozygous mutant. Other symbols will be
#' treated as NA values. A sample raw input file is stored in the
#' \code{inst/extdata} folder of this package. The \code{G2_file} can contain
#' more or less G2 dam mice and more or less genes than the
#' \code{main_file}. The algorithm will calculate the p values for lethality
#' according to different algorithms based on whether G2 dam mice are available
#' for a certain gene. The final p value is the combined p values from all the
#' G3 mice.  A PDF file will be generated. It will contain a diagnostic plot if
#' continuous phenotype scores are provided in the \code{main_file}. It will
#' contain another three manhattan plots, each assuming one of the three
#' "additive", "recessive" and "dominant" models. The two cutoff values are 0.05
#' and 0.05/(number of genes). The p values are transformed to a -log10
#' scale. It will also contain a plot for the analysis result of TDT analysis.
#' A CSV file containing will be generated. Columns names are in the first
#' line. The REF, HET, VAR columns are the number of mice with each genotype for
#' each gene. The "lethal" column is the p value for testing homozygous
#' lethality. The "additive", "recessive" and "dominant" columns are the p
#' values for testing the significance of the correlation of each gene and the
#' phenotype under different models. The "TDT" column is the p value of TDT test
#' for testing each gene's significance. The "Penetrance_REF", "Penetrance_HET",
#' "Penetrance_VAR" columns are the proportions of G3 mice with each genotype
#' showing the "affected" phenotype (only applies to binary phenotype scores or
#' binary scores that are converted from continuous phenotype scores). The
#' "Semidominance" column is the semidominance of the mutations, only calculated
#' when continuous phenotype scores are
#' provided. Semidominance=(HET-REF)/(VAR-REF). 0 means the HET phenotype is the
#' same as or over the average of the REF phenotype and 1 means the HET
#' phenotype is the same as or over the average of the VAR phenotype. Note that
#' VAR may not necessarily be lower than REF.  A second PDF file will be
#' generated containing scattorplots of phenotypes scores (continuous variable)
#' or table plots of phenotypes (binary variable).  A TXT file containing
#' messages, warnings and errors of the program will be generated.  The
#' statistical model used to calculate the p value of lethality, more accurately
#' speaking, is testing "whether homozygous mutation will have any detrimental
#' effect on the survival of mice at the time of data collection". chrX genes
#' are ignored for the moment.
#'
#' @return No return value.
#' @export
#' @examples
#' main_file <-
#'   system.file("extdata/Aquamarine_QM_Report.csv",package="LinkageAnalysis")
#' main_file_continuous <-
#'   system.file("extdata/Aquamarine_QM_Report_continuous.csv",package="LinkageAnalysis")
#' G2_file <-
#'   system.file("extdata/Aquamarine_QM_Report_G2.csv",package="LinkageAnalysis")
#' output <-
#'   sub(pattern="Aquamarine_QM_Report.csv",replacement="output",x=main_file)
#'
#' # if categorical phenotype scores are given
#' single_link(main_file,G2_file,output=output,silent=TRUE,prefix="single")
#' # if continuous numerical phenotype scores are given
#' single_link(main_file_continuous,G2_file,output=output,silent=TRUE,prefix="single")
single_link <- function(vcfFile, pedFile, pheno.name,
                        output = ".", test = "wG2",
                        detect = "never", silent = T, tail = "decreasing",
                        prefix = "", plot.it = TRUE,
                        transform.pheno = NULL) {
  log.file <- filename(output, prefix)$log_file
  ret <- tryCatch(
      {
        ret <- single.link.impl(vcfFile, pedFile, pheno.name,
                                output, test, detect, silent, tail,
                                prefix, plot.it, transform.pheno)
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
        snapshot("single.link.impl", "debug.single.link.impl.Rdata")
        print(str(err))
        print(err)
        msg <- ifelse(is.null(err[["message"]]), "UnknownError", err$message)
        msg <- paste("Exit failed", msg, sep = " ")
        report("m", msg, log.file)
        return(list(returncode = 1, message = msg))
      })

  return(ret)
}

single.link.impl <- function(vcfFile, pedFile, pheno.name,
                             output = ".", test = "wG2",
                             detect = "never", silent = T, tail = "decreasing",
                             prefix = "", plot.it = TRUE,
                             transform.pheno = NULL) {
  start.time <- Sys.time()

  ## set up log file
  log.file <- file.path(getwd(), file.path(output, "log.txt"))
  basicConfig()
  addHandler(writeToFile, file = log.file)
  fns <- filename(output, prefix)  # generate output file names

  # check validity of parameters
  fns <- filename(output, prefix)  # generate output file names
  if (!tail %in% c("increasing", "decreasing", "both")) {
    report("e", "Unrecognized option for tail!", fns$log_file)
  }

  report("m", paste("Version:", packageVersion("LinkageAnalysis")), fns$log_file)
  report("m", paste("Date:", Sys.time()), fns$log_file)
  report("m", paste("Host:", Sys.info()["nodename"]) , fns$log_file)
  report("m", paste("Call:", deparse(sys.status()$sys.calls[[1]])), fns$log_file)


  ## read data
  snapshot("single.link.impl", "debug.single.before.load.Rdata")
  tmp <- load.vcf.ped(vcfFile, pedFile, pheno.name)
  if (isSuccess(tmp)) {
    loginfo("VCF/PED loading succeed.")
    vcf <- tmp$vcf
    ped <- tmp$ped
  } else {
    loginfo("VCF/PED loading failed.")
    return(list(returncode = 1, message = "VCF/PED loading failed."))
  }
  snapshot("single.link", "dbg.single.Rdata")

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
  ret <- data.frame(Gene = gene, chr = vcf$CHROM, pos = vcf$POS,
                    REF = nas ,HET =nas, VAR = nas,
                    lethal = nas, additive = nas, recessive = nas, dominant = nas, TDT = nas,
                    Penetrance_REF = nas, Penetrance_HET = nas, Penetrance_VAR = nas, Semidominance = nas,
                    NMISS = nas, lethalCount = nas)
  rownames(ret) <- NULL

  report("m", "Load data complete", fns$log_file)
  snapshot("single.link.impl", "debug.single.load.Rdata")

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
  null.model <- create.null.model(pheno, pheno.name, test)
  has.random.effect <- grepl("\\(", null.model)
  isBinary <- is.factor(pheno[,pheno.name])

  # fit null model
  null <- fit.null.model(null.model, pheno, isBinary, has.random.effect)
  if (isSuccess(null)) {
    loginfo("Finished fitting null model\n")
  } else {
    logerror("Fitting null model failed!")
    return(null)
  }

  ## start fitting alternative models for each variant
  snapshot("single.link", "dbg.single.fit.alt.Rdata", force = TRUE)

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
      null <- fit.null.model(null.model, pheno, isBinary, has.random.effect)
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
      loginfo("Perform %s test.\n", type)

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
      alt.model <- paste0(null.model, " + gt")
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
        pval <- anova(null, alt)$"Pr(>Chisq)"[2]
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

  snapshot("calc.genetic", "calc.genetic.Rdata")
  ret <- calc.genetic(ret, geno, pheno, pheno.name, isBinary)
  head(ret)

  if (plot.it) {
    snapshot("plot.it", "plot.it.Rdata")
    ## draw linkage plot
    linkage.plot.pdf <- fns$linkage_file
    pdf(file = linkage.plot.pdf, height = 8, width = 20)
    plot.hist(pheno[,pheno.name], main = "Histogram of phenotype scores")
    plot.manhattan(data.frame(Chrom = ret$chr, Position = ret$pos, Gene = ret$Gene, Pval = ret$additive) , main = "Manhattan Plot: additive  model ")
    plot.manhattan(data.frame(Chrom = ret$chr, Position = ret$pos, Gene = ret$Gene, Pval = ret$recessive), main = "Manhattan Plot: recessive model ")
    plot.manhattan(data.frame(Chrom = ret$chr, Position = ret$pos, Gene = ret$Gene, Pval = ret$dominant) , main = "Manhattan Plot: dominant  model ")
    dev.off()
    loginfo(paste0("Generated ", linkage.plot.pdf))

    ## draw distribution plot
    dist.plot.pdf <- fns$distrib_file
    loginfo("%d distribution plots to be generated", length(dist.plots))
    ## require(gridExtra)
    tmp <- do.call(marrangeGrob, c(dist.plots, list(nrow=2, ncol=2)))
    ggsave(dist.plot.pdf, tmp, width = 6, height = 6)
    loginfo(paste0("Generated ", dist.plot.pdf))
  }
  write.table(ret, file = fns$csv_file, quote = F, row.names = F, sep = ",")
  loginfo(paste0("Generated ", fns$csv_file))
  wd <- getwd()
  save(list = ls(), file = fns$results_file)

  end.time <- Sys.time()
  diff.time <- difftime(end.time, start.time, units = "secs")
  msg <- (sprintf("single_link() finished in %.3f seconds - [vcfFile=%s;pedFile=%s;pheno.name=%s;test=%s;detect=%s;tail=%s]",
                  diff.time,
                  vcfFile,
                  pedFile,
                  pheno.name,
                  test,
                  detect,
                  tail))
  print(msg)
  loginfo(msg)
  loginfo("Exit successfully")
  return(list(returncode = 0, message = "", result = ret))
}

## single.link.impl.csv <- function(main_file = "", G2_file = "",
##                              output = ".", test = "wG2",
##                              detect = "never", silent = T, tail = "decreasing",
##                              prefix = "", plot.it = TRUE,
##                              transform.pheno = NULL) {
##   start.time <- Sys.time()

##   fns <- filename(output, prefix)  # generate output file names

##   report("m", paste("Version:", packageVersion("LinkageAnalysis")), fns$log_file)
##   report("m", paste("Date:", Sys.time()), fns$log_file)
##   report("m", paste("Host:", Sys.info()["nodename"]) , fns$log_file)
##   report("m", paste("Call:", deparse(sys.status()$sys.calls[[1]])), fns$log_file)

##   if (!tail %in% c("increasing", "decreasing", "both")) {
##     report("e", "Unrecognized option for tail!", fns$log_file)
##   }

##   # read data
##   snapshot("single.link.impl", "debug.single.before.load.Rdata")
##   tmp <- get_data(main_file, G2_file, fns$log_file, detect, transform.pheno)
##   ## Before data load can be possibily fail, just write pval = 1 as results
##   if (!is.null(tmp$data$genes)) {
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
##   if (!is.null(tmp$data$obs)) {
##   } else {
##     tmp$data$obs <- 0
##   }
##   ## load failes, or no genotype
##   if (tmp$returncode || tmp$data$obs <= 0) {
##     ## no samples are genotypes, so quitting
##     tmp$returncode = 0
##     return(tmp)
##   }

##   raw_data <- tmp$data
##   report("m", "Load data complete", fns$log_file)
##   snapshot("single.link.impl", "debug.single.load.Rdata")

##   genes <- raw_data$genes
##   phenotype <- raw_data$phenotype
##   genotype <- raw_data$genotype
##   bin <- raw_data$bin
##   G2 <- raw_data$G2  # G2 is null is G2 dam genotype data are not available

##   # initialize
##   pt <- phenotype$phenotype
##   sex <- phenotype$sex
##   mother <- phenotype$mother
##   effective <- rep(TRUE, dim(genes)[1])  # whether stat test is carried out on each gene?

##   sig_gene <- as.data.frame(matrix(data = 1, ncol = 4, nrow = dim(genes)[1]))  # significance of genes are stored in this vector
##   colnames(sig_gene) <- c("additive", "recessive", "dominant", "TDT")

##   if (plot.it) {
##     # statistical testing and draw the distribution plot at the same time
##     pdf(file = fns$distrib_file, width = 6, height = 9)
##     par(mfrow = c(3, 2), mar = c(4, 4, 4, 4))
##   }

##   for (i in 1:dim(genes)[1]) {
##     # get genotype in character string
##     if (silent == F) {
##       msg <- paste("---", genes$Gene[i], "[", i, "of", dim(genes)[1], "]", "---")
##       report("m", msg, fns$log_file)
##     }
##     gt <- unlist(genotype[i, ])
##     if (sum(!is.na(convert_gt(gt, "additive"))) < 10) {
##       if (silent == F) {
##         report("m", "Skip because there are too few valid values", fns$log_file)
##       }
##       effective[i] <- FALSE
##       next
##     }
##     if (silent == F) {
##       report("m", paste("# remaining mice", sum(!is.na(convert_gt(gt, "additive")))),
##              fns$log_file)
##     }

##     # glmer anova test of full model and reduced model for each of the three types
##     for (type in c("additive", "recessive", "dominant")) {
##       data <- data.frame(pt = pt, sex = sex, gt = convert_gt(gt, type), mother = mother)
##       data <- data[!is.na(data$gt), ]
##       if (length(unique(data$gt)) > 1) {
##         ## use tryCatch to avoid crashing
##         pval <- tryCatch(
##             {
##               pval <- anova_test(data, bin, test, silent, fns$log_file, tail)$pvalue
##             },
##             error = function(err) {
##               snapshot("single.link.impl", "debug.single.link.impl.Rdata")
##               print(str(err))
##               print(err)
##               msg <- ifelse(is.null(err[["message"]]), "UnknownError", err$message)
##               msg <- paste("Fitting failed", msg, sep = " ")
##               report("m", msg, fns$log_file)
##               return(list(returncode = 1, message = msg))
##             })
##         if (is.list(pval) && pval$returncode == 1) {
##           sig_gene[i, type] <- pval
##         } else {
##           sig_gene[i, type] <- pval
##         }
##         if (sig_gene[i,type] == 0) {
##           sig_gene[i,type]=1e-300
##         }
##       }
##     }

##     # TDT analysis
##     data <- data.frame(pt = pt, sex = sex, gt = convert_gt(gt, "additive"), mother = mother)
##     data <- data[!is.na(data$gt), ]
##     sig_gene[i, "TDT"] <- TDT(data, G2, genes[i, ], bin, fns$log_file)

##     # distribution plot
##     if (length(unique(data$gt)) > 1 && plot.it) {
##       distrib_plot(data, genes[i, ], bin)
##     }
##   }

##   if (plot.it) {
##     dev.off()  #close the distribution plot
##     cat("Distribution plot saved at: ", fns$distrib_file, "\n")
##     ## par(mfrow = c(1, 1))
##   }

##   # flag lethal genes
##   genes$REF <- apply(genotype, 1, function(x) sum(x == "REF"))
##   genes$HET <- apply(genotype, 1, function(x) sum(x == "HET"))
##   genes$VAR <- apply(genotype, 1, function(x) sum(x == "VAR"))

##   if (silent == F) {
##     report("m", "Flag lethal genes\n", fns$log_file)
##   }
##   genes$lethal <- single_lethal(genotype, genes, phenotype, G2)

##   alpha <- 0.05
##   bonferroni <- alpha / sum(effective)
##   if (plot.it) {
##     # draw manhattan plot and diagnostic plot
##     if (silent == F) {
##       report("m", "Plotting Manhattan plots\n", fns$log_file)
##     }
##     pdf(file = fns$pdf_file, height = 8, width = 16)
##     par(mfrow = c(1, 1))

##     diagnostic(phenotype, bin, fns$log_file, raw_data$unconverted)  # draw diagnostic plot
##     # manhattan plot
##     for (type in c("additive", "recessive", "dominant", "TDT")) {
##       print(head(genes))
##       ## browser()
##       mht(genes, sig_gene, type, fns$log_file, test, bin, effective,
##           genotype)
##     }
##     dev.off()
##   }
##   # calculate penetrance and semidominance
##   genetics <- get_genetics(genes, phenotype, genotype, bin, raw_data$unconverted)

##   # save results
##   if (silent == F) {
##     msg <- sprintf("Save results in CSV format - %s", fns$log_file)
##     report("m", msg, fns$log_file)
##   }
##   results <- cbind(genes[, c("Gene", "chr", "pos", "REF", "HET", "VAR", "lethal")],
##                    sig_gene, genetics)
##   write.table(results, file = fns$csv_file, quote = F, row.names = F, sep = ",")

##   analysis <- list(test = test, detect = detect, tail = tail, prefix = prefix,
##                    bin = bin, results = results, bonferroni = bonferroni)
##   save(analysis, file = fns$results_file)

##   end.time <- Sys.time()
##   diff.time <- difftime(end.time, start.time, units = "secs")
##   msg <- (sprintf("single_link() finished in %.3f seconds - [main=%s;G2=%s;test=%s;detect=%s;tail=%s]",
##                   diff.time,
##                   main_file,
##                   G2_file,
##                   test,
##                   detect,
##                   tail))
##   print(msg)
##   report("m", msg, fns$log_file)

##   return(list(returncode = 0, result = results))
## }

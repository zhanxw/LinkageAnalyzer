#' @title Statistical analysis for correlation with phenotype and lethality on a single locus
#'
#' @description
#' This function reads the raw data and conducts statistical
#' analysis for the significance of correlation with phenotype and homozygous
#' lethality for each gene.
#'
#' @param vcfFile genotype input file in VCF format
#' @param pedFile pedigree file in PED format
#' @param pheno.name phenotype name to analyze
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
#' @param log.level Default WARN, but can be set from 'DEBUG', 'INFO', 'WARN'
#' and 'ERROR'.
#'
#' @details single.link takes a VCF file as genotype covariates and a PED format
#' as phenotype (response).  A sample input file is stored in the
#' \code{inst/extdata/single} folder of this package. The algorithm will
#' calculate the p values for lethality conditioning on G2 dam mice genotype
#' for a certain gene. The final p value is the combined p values from all the
#' G3 mice.  A PDF file will be generated. It will contain a diagnostic plot if
#' continuous phenotype scores are provided. It will then contain three
#' manhattan plots, each assuming one of the three "additive", "recessive" and
#' "dominant" models. The two cutoff values are 0.05 and 0.05/(number of genes).
#' The p values are transformed to a -log10 scale for better visualization.
#' It also contain a plot for the analysis results from TDT analysis.
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
#' @return a list of two values will be returned. One is returncode (0: success).
#' The other is analysis results.
#' @export
#' @examples
#' path <- system.file("extdata/single",package="LinkageAnalysis")
#' vcfFile <- file.path(path, "R0491_body_weight.vcf")
#' pedFile <- file.path(path, "R0491_body_weight.ped")
#' pheno.name <- "weight"
#' output <- file.path(path, "output")
#' ret <- single.link(vcfFile, pedFile, pheno.name, output, test = "woG2", tail = "both", prefix="single")
single.link <- function(vcfFile,
                        pedFile,
                        pheno.name,
                        output = ".",
                        test = "wG2",
                        detect = "never",
                        silent = T,
                        tail = "decreasing",
                        prefix = "",
                        plot.it = TRUE,
                        transform.pheno = NULL,
                        log.level = 'WARN') {
  collectUsage("single.link")
  log.file <- filename(output, prefix)$log_file
  ret <- tryCatch(
      {
        ret <- single.link.impl(vcfFile, pedFile, pheno.name,
                                output, test, detect, silent, tail,
                                prefix, plot.it, transform.pheno, log.level)
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
        reportError(err)
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
                             transform.pheno = NULL,
                             log.level = 'WARN') {
  start.time <- Sys.time()

  ## set up log file
  log.file <- file.path(getwd(), file.path(output, "log.txt"))
  basicConfig(log.level)
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

  loginfo("Load data complete")
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
  null.model <- create.null.model(pheno, pheno.name, test)
  if (!isSuccess(null.model)) {
    return(null.model)
  }
  has.random.effect <- grepl("\\(", null.model)
  isBinary <- is.factor(pheno[,pheno.name])

  snapshot("calc.genetic", "calc.genetic.Rdata")
  ret <- calc.genetic(ret, geno, pheno, pheno.name, isBinary)
  head(ret)

  # fit null model
  null <- fit.null.model(null.model, pheno, isBinary, has.random.effect)
  if (isSuccess(null)) {
    loginfo("Finished fitting null model\n")
  } else {
    logerror("Fitting null model failed!")
    write.table(ret, file = fns$csv_file, quote = F, row.names = F, sep = ",")
    loginfo(paste0("Generated %s anyway", fns$csv_file))
    return(null)
  }

  ## start fitting alternative models for each variant
  snapshot("single.link", "dbg.single.fit.alt.Rdata")

  dist.data <- list()
  dist.plots <- list()
  for (i in 1:nVariant) {
    if (i > 10 && is.debug.mode()) {
      cat("DEBUG skipped ", i, "th variant ..\n")
      next
    }
    loginfo("Process %d out of %d variant: %s", i, nVariant, gene[i])

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
    save.dist.plot(dist.plots, nrow = 2, ncol = 2, dist.plot.pdf)
    loginfo(paste0("Generated ", dist.plot.pdf))
  }
  write.table(ret, file = fns$csv_file, quote = F, row.names = F, sep = ",")
  loginfo(paste0("Generated ", fns$csv_file))
  wd <- getwd()
  save(list = ls(), file = fns$results_file)

  end.time <- Sys.time()
  diff.time <- difftime(end.time, start.time, units = "secs")
  msg <- (sprintf("single.link() finished in %.3f seconds - [vcfFile=%s;pedFile=%s;pheno.name=%s;test=%s;detect=%s;tail=%s]",
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


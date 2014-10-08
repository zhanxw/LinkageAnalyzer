#' Meta-analyze single variant in super pedigree
#'
#' Perform variant-base superpedigree analysis
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
#' path <- system.file("extdata/meta",package="LinkageAnalysis")
#' vcfFile <- file.path(path, "4_130378043_C_G.mcmv.vcf")
#' pedFile <- file.path(path, "4_130378043_C_G.mcmv.ped")
#' pheno.name <- "log"
#' output <- file.path(path, "output")
#' ret <- meta.single.link(vcfFile, pedFile, pheno.name, output, test = "woG2", tail = "both", prefix="meta")
meta.single.link <- function(vcfFile, ## a vector of list
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
  collectUsage("meta.single.link")
  log.file <- filename(output, prefix)$log_file
  ret <- tryCatch(
      {
        ret <- meta.single.link.impl(vcfFile, pedFile, pheno.name,
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
        snapshot("meta.single.link.impl", "debug.meta.single.link.impl.Rdata")

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

run.mixed.effect.alt.model <- function(isBinary, alt.model, pheno) {
  assign("last.warning", NULL, envir = baseenv())
  alt <- tryCatch(
      {
        if (isBinary) {
          alt <- glmer(as.formula(alt.model), data = pheno, family = "binomial")
        } else {
          alt <- lmer(as.formula(alt.model), data = pheno, REML = FALSE)
        }
      },
      warning = function(warn) {
        logwarn(paste0("Fit [ ", alt.model, " binary =", isBinary, " ] has warnings ", str(warn)))
        msg <- as.character(get("last.warning", baseenv()))
        return(list(returncode = 1, message = msg, warning = warn, isBinary = isBinary, alt.model = alt.model))
      },
      error = function(err) {
        logwarn(paste0("Fit [ ", alt.model, " binary =", isBinary, " ] failed ", str(err)))
        reportError(err)
        msg <- ifelse(is.null(err[["message"]]), "UnknownError", err$message)
        return(list(returncode = 1, message = msg, error = err, isBinary = isBinary, alt.model = alt.model))
      })
  alt
}

run.fixed.effect.alt.model <- function(isBinary, alt.model, pheno, pheno.name, tail) {
  fid <- numGT <- NULL ## bypass CRAN check
  tmp <- ddply(pheno, .(fid), function(x) {
    c(numGT = length(unique(x$gt)))})
  tmp <- subset(tmp, numGT == 1)$fid
  reduced.pheno <- subset (pheno, !fid %in% tmp )
  reduced.model <- str_replace(alt.model, "\\(1\\|fid\\)", "fid")
  if (length(unique(reduced.pheno$fid)) == 1) {
    reduced.model <- str_replace(reduced.model, "\\+ fid", "")
  }
  reduced.model <- str_replace(reduced.model, "\\(1\\|mother\\)", "mother")
  if (length(unique(reduced.pheno$mother)) == 1) {
    reduced.model <- str_replace(reduced.model, "\\+ mother", "")
  }
  reduced.pheno$mother <- factor(reduced.pheno$mother)

  ## check if the model can be fit
  model.fittable <- TRUE
  while (model.fittable) {
    ## check sample size
    if (nrow(reduced.pheno) == 0) {
      loginfo("Insufficient sample size to fit models, set pvalue as one")
      pval <- 1
      model.fittable <- FALSE
      break
    }

    ## check factor variables
    model.var <- str_split(reduced.model, " ")[[1]]
    model.var <- model.var [!model.var %in% c("", "~", "1", "+") ]
    model.var <- reduced.pheno[, model.var]
    model.var <- model.var[complete.cases(model.var), ]
    useless.var <- apply(model.var, 2, function(x) {length(unique(x)) == 1})

    if (any(useless.var)) {
      loginfo(paste0("Reduced model has useless variable:", names(model.var)[useless.var], ", set pvalue to one."))
      pval <- 1
      model.fittable <- FALSE
      break
    }

    ## check rank
    reduced.model.matrix <- model.matrix(as.formula(reduced.model), reduced.pheno)
    if ( qr(reduced.model.matrix)$rank < ncol(reduced.model.matrix)) {
      loginfo(paste0("Insufficient sample size to fit models, set pvalue to one."))
      pval <- 1
      model.fittable <- FALSE
      break
    }

    ## try alt model
    if (isBinary) {
      reduced.pheno[, pheno.name] <- as.numeric(reduced.pheno[, pheno.name]) - 1
      alt <- glm(as.formula(reduced.model), data = reduced.pheno, family = "binomial")
    } else {
      alt <- lm(as.formula(reduced.model), data = reduced.pheno)
    }
    idx <- which(colnames(summary(alt)$coefficients) == "Pr(>|t|)" |
                 colnames(summary(alt)$coefficients) == "Pr(>|z|)")
    coef <- summary(alt)$coefficients
    if ("gt" %in% rownames(coef) ) {
      pval <- coef["gt", idx]
      direction <- coef(alt)["gt"] > 0
      if (!is.na(direction)) {
        pval <- convert_tail(direction, pval, tail)
      }
    } else{
      ## glm fitting failed, so set pval to one
      logwarn("Refit using glm() failed, set pvalue to one.\n")
      pval <- 1
    }
    break
  }
  pval
}

create.null.model <- function(pheno, pheno.name, test) {
  null.model <- sprintf("%s ~ 1 ", pheno.name)
  if (length(unique(pheno$fid)) > 1) {
    null.model <- paste0(null.model, " + (1|fid) ")
  } else {
    loginfo("INFO: 1 family detected\n")
  }
  if (length(unique(pheno$sex)) > 1) {
    null.model <- paste0(null.model, " + sex")
  } else {
    loginfo("INFO: only one gender detected\n")
  }

  if (test == "wG2") {
    if (all(is.na(pheno$mother))) {
      msg <- "Does not have mother info at all (wG2 cannot work)!!"
      logerror(msg)
      return(list(returncode = 1, message = msg))
    }
    if (length(unique(pheno$mother))==1) {
      msg <- "Cannot model G2 effect as there is only one G2 mouse."
      logerror(msg)
      return(list(returncode = 1, message = msg))
    }
    null.model <- paste0(null.model, " + (1|mother)")
  }
  null.model
}

isResponseContant <- function(model, pheno) {
  pheno.name = str_extract(model, "^[^~ ]+")
  y <- pheno[, pheno.name]
  if (length(y) == 0 || length(unique(y)) == 1) {
    return(TRUE)
  }
  return(FALSE)
}

countModelParam <- function(model, pheno, isBinary, has.random.effect) {
  if (has.random.effect) {
    nobars(model)
    num.fixed.eff <- ncol(model.matrix(nobars(as.formula(model)), data = pheno))
    num.random.eff <- do.call(sum, lapply(findbars(as.formula(model)),
                             function(x) {
                               eff <- (sub("1 \\| ", "", deparse(x)))
                               length(unique(pheno[, eff]))
                             }))
    return(num.fixed.eff + num.random.eff)
  }

  num.fixed.eff <- ncol(model.matrix(as.formula(model), data = pheno))
  return(num.fixed.eff)
}

fit.null.model <- function(null.model, pheno, isBinary, has.random.effect) {
  loginfo("Fit null model: %s", null.model)
  snapshot("fit.null.model", "fit.null.model.Rdata")
  if (isResponseContant(null.model, pheno)) {
    return(list(returncode = 1, message = "Response is constant - cannot fit the model"))
  }
  if (grepl("mother", null.model) && length(unique(pheno$mother)) == 1)  {
    return(list(returncode = 1, message = "Cannot access mother effect when input only contains one mother mouse"))
  }
  nSample <- nrow(pheno)
  nParam <- countModelParam(null.model, pheno, isBinary, has.random.effect)
  if (nSample <= nParam) {
    return(list(returncode = 1, message = "Sample size is smaller than free model parameters"))
  }
  if (isBinary) {
    loginfo("Phenotypes are treated as binary\n")
    null <- tryCatch(
        {
          if (has.random.effect) {
            glmer(as.formula(null.model), data = pheno, family = "binomial")
          } else {
            glm(as.formula(null.model), data = pheno, family = "binomial")
          }
        },
        error = function(err) {
          reportError(err)
          msg <- ifelse(is.null(err[["message"]]), "UnknownError", err$message)
          return(list(returncode = 1, message = msg, error = err))
        })
    if (!isSuccess(null)) {
      null.model <- str_replace(null.model, "\\+ sex", "")
      loginfo("Refit null model using less covariates: %s\n", null.model)
      null <- tryCatch(
          {
            if (has.random.effect) {
              null <- glmer(as.formula(null.model), data = pheno, family = "binomial")
            } else {
              null <- glm(as.formula(null.model), data = pheno, family = "binomial")
            }
            null
          },
          error = function(err) {
            reportError(err)
            msg <- ifelse(is.null(err[["message"]]), "UnknownError", err$message)
            return(list(returncode = 1, message = msg, error = err))
          })
    }
  } else { ## qtl
    loginfo("Phenotypes are treated as continuous\n")
    null <- tryCatch(
        {
          if (has.random.effect) {
            null <- lmer(as.formula(null.model), data = pheno, REML = FALSE)
          } else {
            null <- lm(as.formula(null.model), data = pheno)
          }
          null
        },
        error = function(err) {
          reportError(err)
          msg <- ifelse(is.null(err[["message"]]), "UnknownError", err$message)
          return(list(returncode = 1, message = msg, error = err))
        })
  }
  null
}

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

meta.single.link.impl <- function(vcfFile, ## a vector of list
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
  basicConfig(log.level)
  addHandler(writeToFile, file = log.file)
  ## if (file.exists(log.file)) {
  ##   file.remove(log.file)
  ## }
  fns <- filename(output, prefix)  # generate output file names

  ## examine input parameters
  if (!tail %in% c("increasing", "decreasing", "both")) {
    msg <- paste0("ERROR: Unrecognized option for tail [ ", tail, "]\n")
    logerror(msg)
    return(list(returncode = 1, message = msg))
  }

  ## record running environment
  wd <- getwd()
  loginfo(paste("Version:", packageVersion("LinkageAnalysis")))
  loginfo(paste("Date:", Sys.time()))
  loginfo(paste("Host:", Sys.info()["nodename"]))
  loginfo(paste("Call:", deparse(sys.status()$sys.calls[[1]], width.cutoff = 500L)))
  loginfo(paste("Directory:", wd))


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
  snapshot("meta.link", "dbg.meta.Rdata")

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
    loginfo("Generated %s anyway", fns$csv_file)
    return(null)
  }

  ## start fitting alternative models for each variant
  snapshot("meta.link", "dbg.meta.fit.alt.Rdata")

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
    save.dist.plot(dist.plots, nrow = 2, ncol = 1, fn = dist.plot.pdf)
    loginfo(paste0("Generated ", dist.plot.pdf))
  }
  write.table(ret, file = fns$csv_file, quote = F, row.names = F, sep = ",")
  loginfo(paste0("Generated ", fns$csv_file))
  wd <- getwd()
  save(list = ls(), file = fns$results_file)

  end.time <- Sys.time()
  diff.time <- difftime(end.time, start.time, units = "secs")
  msg <- (sprintf("meta.single.link() finished in %.3f seconds - [vcfFile=%s;pedFile=%s;pheno.name=%s;test=%s;detect=%s;tail=%s]",
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

check.mixed.effect.alt.model.converge <- function(alt) {
  if (.hasSlot(alt, "optinfo")) {
    converge =  is.null(alt@optinfo$conv$lme4$messages[[1]])
    if (!converge) {
      msg <- alt@optinfo$conv$lme4$messages[[1]]
      err <- list(message = msg)
      class(err) <- "error"
      logwarn(paste0("Model may suffer from convergence: ", msg))
      alt <- list(returncode = 1, message = msg, error = err)
      return(alt)
    }
  }
  return(alt)
}

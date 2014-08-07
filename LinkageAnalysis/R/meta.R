#' Meta-analyze single variant in super pedigree
#'
#' Perform single analysis using multiple pedigrees
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
#'
#' @export
#' @examples
#' vcfFile <-
#'   system.file("extdata/R0359-R0360.vcf",package="LinkageAnalysis")
#' pedFile <-
#'   system.file("extdata/R0359-R0360.ped",package="LinkageAnalysis")
#' pheno.name <- "raw"
#' output <-
#'   sub(pattern="R0359-R0360.vcf",replacement="output",x=vcfFile)
#' meta.single.link(vcfFile, pedFile, pheno.name, output = output)
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
                             transform.pheno = NULL) {
  log.file <- filename(output, prefix)$log_file
  ret <- tryCatch(
      {
        ret <- meta.single.link.impl(vcfFile, pedFile, pheno.name,
                                     output, test,
                                     detect, silent, tail, prefix,
                                     plot.it, transform.pheno)
        if (ret$returncode == 0) {
          msg <- paste("Exit successfully", ret$message, sep = " ")
        } else {
          ## deal with ignorable errors
          if (ret$returncode == 1 && ret$message == "dichotomize failed") {
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

        print(str(err))
        print(err)
        msg <- ifelse(is.null(err[["message"]]), "UnknownError", err$message)
        msg <- paste("Exit failed", msg, sep = " ")
        report("m", msg, log.file)
        return(list(returncode = 1, message = msg))
      })

  return(ret)
}

#' Examine mother-offspring pair and fix mother's genotype by her offsprings.
#' e.g. mom has VAR offsprings, mom genotype cannot be REF => mom must be HET
#' @param pheno ped data
crossCheckMotherGenotype <- function(pheno) {
  mother <- NULL # bypass CRAN warning
  mom <- ddply(pheno, .(mother), function(x) {
    c(
        ref = sum(x$gt == 0, na.rm = TRUE),
        var = sum(x$gt > 0, na.rm = TRUE))
  })
  total.fix <- 0
  fixed.mom.id <- vector("character", 0)
  for (i in seq_len(nrow(mom))) {
    mom.name <- mom[i, "mother"]
    if (mom[i, "var"] > 0) {
      if (!mom.name %in% pheno$iid ) {
        next
      }
      if (!is.na(pheno[ pheno$iid == mom.name, "gt"]) &&
          pheno[ pheno$iid == mom.name, "gt"] != 1) {
        pheno[ pheno$iid == mom.name, "gt"] <- 1 ## set to het
        total.fix <- total.fix + 1
        fixed.mom.id <- c(fixed.mom.id, mom.name)
      }
    }
  }
  if (total.fix) {
    cat("INFO: Fixed ", total.fix, " mother genotypes: ", unique(fixed.mom.id), "\n")
  }
  pheno
}

countForLethal <- function(pheno) {
  nHetMom <- nVarFromHetMom <- nUnknownMom <- nVarFromUnknownMom <- 0
  mothers <- unique(pheno[, "mother"])
  g2 <- pheno$iid[pheno$gen == 2]
  g2.mothers <- intersect(g2, mothers)
  if (length(mothers) == 0) {
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
      print(mother)
      print(offspring)
    } else if (mother.gt == 0) {
      next
    } else if (mother.gt == 1) {
      offspring <- pheno[pheno$mother == mother & !is.na(pheno$gt), ]
      if (nrow(offspring) > 0 ) {
        nHetMom        <- nHetMom        + sum(!is.na(offspring$gt), na.rm = TRUE)
        nVarFromHetMom <- nVarFromHetMom + sum(offspring$gt == 2, na.rm = TRUE)
      }
      print(mother)
      print(offspring)
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
                                  transform.pheno = NULL) {
  ## create log file
  log.file <- file.path(output, "log.txt")
  if (file.exists(log.file)) {
    file.remove(log.file)
  }

  mycat <- function(...) {
    cat(...)
    cat(..., file = log.file, append = TRUE)
  }

  wd <- getwd()
  mycat(paste("Version:", packageVersion("LinkageAnalysis")), "\n")
  mycat(paste("Date:", Sys.time()), "\n")
  mycat(paste("Host:", Sys.info()["nodename"]), "\n")
  mycat(paste("Call:", deparse(sys.status()$sys.calls[[1]], width.cutoff = 500L)), "\n")
  mycat(paste("Directory:", wd), "\n")

  start.time <- Sys.time()
  mycat("Load VCF file: ", vcfFile, "\n")
  vcf <- get.vcf(vcfFile)
  vcf.summarize(vcf)

  mycat("Load PED file: ", pedFile, "\n")
  ped <- get.ped(pedFile)
  ped.summarize(ped)

  ## verify phenotype name
  stopifnot(pheno.name %in% colnames(ped)[-(1:5)] )
  stopifnot(!all(is.na(ped[,pheno.name])))

  ## clean up samples in VCF and PED
  idx <- ! vcf$sampleId %in% ped[,2]
  mycat("Remove ", sum(idx), " samples from VCF as they are not in PED.\n")  ## some sample may not be screened
  vcf <- vcf.delete.sample.by.index(vcf, idx)

  tmp <- setdiff(ped$iid, vcf$sampleId)
  mycat("Add ", length(tmp), " samples to VCF according to PED file.\n")
  vcf <- vcf.add.sample(vcf, tmp)

  ## rearrange PED according to VCF
  stopifnot(all(sort(vcf$sampleId) == sort(ped$iid)))
  idx <- match(vcf$sampleId, ped$iid)
  ped <- ped[idx,]
  stopifnot(all(vcf$sampleId == ped$iid))

  mycat("VCF/PED loaded\n")
  snapshot("meta.link", "dbg.meta.Rdata")

  mycat("Detect phenotype type\n")
  if (detect == "auto") {
    tmp <- dichotomize(ped[,pheno.name])
    if (tmp$succ) {
      ped[,pheno.name] <- tmp$new.value
    } else {
      mycat("ERROR: Dichotomize failed\n")
      # write status.file
      status.file.name <- file.path(dirname(log.file),
                                    "R_jobs_complete_with_no_output.txt")
      cat(date(), file = status.file.name)
      cat("\t", file = status.file.name, append = TRUE)
      # write log
      msg <- sprintf("Log file [ %s ] created.", status.file.name)
      mycat(msg)
      msg <- "dichotomize failed"
      mycat(msg)
      return(list(returncode = 1, message = msg))
    }
  } else {
    ## phenotype is quantitative, but let's double check it's QTL
    tmp <- ped[,pheno.name]
    tmp <- tmp[!is.na(tmp)] ## remove NA
    tmp <- tmp[! tmp %in% c(-9, 0, 1, 2) ]
    if (length(unique(tmp)) == 0) {
      mycat("WARNING: phenotype seems bo be binary and will apply binary models\n")
      tmp <- ped[,pheno.name]
      idx <- tmp %in% c(-9, 0)
      ped[idx, pheno.name] <- NA
      ped[,pheno.name] <- factor(ped[,pheno.name], levels = c(1, 2), labels = c("AFFECTED", "UNAFFECTED"))
    }
  }

  fns <- filename(output, prefix)  # generate output file names
  if (!tail %in% c("increasing", "decreasing", "both")) {
    msg <- paste0("ERROR: Unrecognized option for tail [ ", tail, "]\n")
    mycat(msg)
    return(list(returncode = 1, message = msg))
  }

  ## make an output skeleton
  nSite <- length(vcf$CHROM)
  nas <- rep(NA, nSite)
  gene <- str_replace(sapply(vcf$INFO, function(x) {str_split(x, ";")[[1]][1]}), "GENE=", "")
  ret <- data.frame(Gene = gene, chr = vcf$CHROM, pos = vcf$POS,
                    REF = nas ,HET =nas, VAR = nas,
                    lethal = nas, additive = nas, recessive = nas, dominant = nas, TDT = nas,
                    Penetrance_REF = nas, Penetrance_HET = nas, Penetrance_VAR = nas, Semidominance = nas)
  rownames(ret) <- NULL

  # calculate lethal
  type <- "lethal"
  mycat("Perform ", type, " test.\n")
  for (i in seq_len(nrow(ret))) {
    ## calculate lethal
    stopifnot(all(vcf$sampleId == ped$iid))
    ## encode genotypes
    tmp <- ped
    tmp$gt <- convert_gt(vcf$GT[i,], "additive")
    ## check mother genotypes
    tmp <- crossCheckMotherGenotype(tmp)
    ## count mother, offspring by their genotypes
    tmp <- countForLethal(tmp)
    ret[i, "lethal"] <- do.call(single.lethal.getPvalue, tmp)
  }

  # reduce data by:
  #  1. select G3 mice
  #  2. remove mice with missing phenotype
  # get G3 mice
  pheno <- ped[!is.na(ped[,pheno.name]), ]
  pheno <- pheno[!is.na(pheno$gen), ]
  pheno <- pheno[pheno$gen == 3, ]

  # prepare genotype data
  geno <- vcf$GT[, pheno$iid]  # G3 genotypes
  nVariant <- nrow(geno)

  null.model <- sprintf("%s ~ 1 ", pheno.name)
  if (length(unique(pheno$fid)) > 1) {
    null.model <- paste0(null.model, " + (1|fid) ")
  } else {
    mycat("INFO: 1 family detected\n")
  }
  if (length(unique(pheno$sex)) > 1) {
    null.model <- paste0(null.model, " + sex")
  } else {
    mycat("INFO: only one gender detected\n")
  }

  if (grepl("\\(", null.model)) {
    has.random.effect <- TRUE
  } else {
    has.random.effect <- FALSE
  }

  if (test == "wG2") {
    if (all(is.na(pheno$mother))) {
      mycat("Does not have mother info at all (wG2 cannot work)!!\n")
      stop("Quit...")
    }
    null.model <- paste0(null.model, " + (1|mother)")
  }

  isBinary <- is.factor(pheno[,pheno.name])
  ## require(lme4)
  if (isBinary) {
    mycat("INFO: Phenotypes are treated as binary\n")
    null <- tryCatch(
        {
          if (has.random.effect) {
            glmer(as.formula(null.model), data = pheno, family = "binomial")
          } else {
            glm(as.formula(null.model), data = pheno, family = "binomial")
          }
        },
        error = function(err) {
          print(str(err))
          print(err)
          msg <- ifelse(is.null(err[["message"]]), "UnknownError", err$message)
          return(list(returncode = 1, message = msg, error = err))
        })
    if (is.list(null) && null$returncode == 1) {
      mycat("Refit null model using less covariates\n")
      null.model <- str_replace(null.model, "\\+ sex", "")
      if (has.random.effect) {
        null <- glmer(as.formula(null.model), data = pheno, family = "binomial")
      } else {
        null <- glm(as.formula(null.model), data = pheno, family = "binomial")
      }
    }
  } else { ## qtl
    mycat("INFO: Phenotypes are treated as continuous\n")
    if (has.random.effect) {
      null <- lmer(as.formula(null.model), data = pheno, REML = FALSE)
    } else {
      null <- lm(as.formula(null.model), data = pheno)
    }
  }
  mycat("Finished fitting null model\n")

  snapshot("meta.link", "dbg.meta.fit.alt.Rdata")
  dist.plots <- list()
  for (i in 1:nVariant) {
    ## for (i in 1:5) {
    if (i > 10 && is.debug.mode()) {
      cat("DEBUG skipped ", i, "th variant ..\n")
      next
    }
    mycat("Process ", i, " out of ", nVariant, " variant : ", gene[i], "\n")
    if (length(unique(geno[i,])) == 1) {
      ## skip mono site
      next
    }
    for (type in c("additive", "recessive", "dominant")) {
      mycat("Perform ", type, " test.\n")
      ## encode genotypes
      pheno$gt <- convert_gt(geno[i,], type)

      ## impute missing genotypes
      pheno$gt[is.na(pheno$gt)] <- 2

      ## skip this site when the genotypes are monomorphic
      if (length(unique(pheno$gt)[!is.na(unique(pheno$gt))]) <= 1) {
        mycat("Skip monomorphic site under ", type, " model\n")
        ret[i, type] <- pval <- 1
        next
      }

      # record graph
      if (type == "additive") {
        ## cat ("store graph\n")
        dist.title <- sprintf("%s (%s:%s)", gene[i], as.character(ret$chr[i]), as.character(ret$pos[i]))
        dist.plots[[length(dist.plots) + 1]] <- plot.distribution(pheno, dist.title)
      }

      ## fit alternative model
      alt.model <- paste0(null.model, " + gt")
      run.mixed.effect.alt.model <- function(isBinary, alt.model, pheno) {
        alt <- tryCatch(
            {
              if (isBinary) {
                alt <- glmer(as.formula(alt.model), data = pheno, family = "binomial")
              } else {
                alt <- lmer(as.formula(alt.model), data = pheno, REML = FALSE)
              }
            },
            warning = function(warn) {
              mycat("Fit [ ", alt.model, " binary=", isBinary, " ] has warnings ", str(err), "\n")
              msg <- as.character(last.warning)
              return(list(returncode = 1, message = msg, warning = warn, isBinary = isBinary, alt.model = alt.model))
            },
            error = function(err) {
              mycat("Fit [ ", alt.model, " binary=", isBinary, " ] failed ", str(err), "\n")
              ## print(err)
              msg <- ifelse(is.null(err[["message"]]), "UnknownError", err$message)
              return(list(returncode = 1, message = msg, error = err, isBinary = isBinary, alt.model = alt.model))
            })
        alt
      }
      assign("last.warning", NULL, envir = baseenv())
      if (has.random.effect) {
        alt <- run.mixed.effect.alt.model(isBinary, alt.model, pheno)
      } else {
        alt <- list(returncode = 1, message = "no random effect")
      }

      ## check convergence
      if (.hasSlot(alt, "optinfo")) {
        converge =  is.null(alt@optinfo$conv$lme4$messages[[1]])
        if (!converge) {
          msg <- alt@optinfo$conv$lme4$messages[[1]]
          err <- list(message = msg)
          class(err) <- "error"
          mycat("Model may suffer from convergence: ", msg, "\n")
          alt <- list(returncode = 1, message = msg, error = err, isBinary = isBinary, alt.model = alt.model)
        }
      }

      ## calculate p-value
      if (!(is.list(alt) && alt$returncode == 1)) {
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
        mycat("Refit alternative model using reduced data and Wald test\n")

        run.fixed.effect.alt.model <- function(isBinary, alt.model, pheno) {
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
              mycat("Insufficient sample size to fit models, set pvalue as one\n")
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
              mycat("Reduced model has unless variable:", names(model.var)[useless.var], ", set pvalue to one.\n")
              pval <- 1
              model.fittable <- FALSE
              break
            }

            ## check rank
            reduced.model.matrix <- model.matrix(as.formula(reduced.model), reduced.pheno)
            if ( qr(reduced.model.matrix)$rank < ncol(reduced.model.matrix)) {
              mycat("Insufficient sample size to fit models, set pvalue to one.\n")
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
              mycat("Refit using glm() failed, set pvalue to one.\n")
              pval <- 1
            }
            break
          }
          pval
        }
        pval <- run.fixed.effect.alt.model(isBinary, alt.model, pheno)

      }
      ret[i, type] <- pval
    }
    ret[i, "TDT"] <- NA     ## TODO: implement this
  }

  ret$REF <- rowSums(geno == 0, na.rm = TRUE)
  ret$HET <- rowSums(geno == 1, na.rm = TRUE)
  ret$VAR <- rowSums(geno == 2, na.rm = TRUE)

  if (isBinary) {
    ret$Penetrance_REF <- apply(geno, 1, function(x) {idx <- x==0; mean(pheno[idx, pheno.name] == "AFFECTED", na.rm = TRUE)})
    ret$Penetrance_HET <- apply(geno, 1, function(x) {idx <- x==1; mean(pheno[idx, pheno.name] == "AFFECTED", na.rm = TRUE)})
    ret$Penetrance_ALT <- apply(geno, 1, function(x) {idx <- x==2; mean(pheno[idx, pheno.name] == "AFFECTED", na.rm = TRUE)})
  }
  ret$Semidominance <- apply(geno, 1, function(x) {
    idx <- x == 0; m0 <- mean(as.numeric(pheno[idx, pheno.name]), na.rm = TRUE)
    idx <- x == 1; m1 <- mean(as.numeric(pheno[idx, pheno.name]), na.rm = TRUE)
    idx <- x == 2; m2 <- mean(as.numeric(pheno[idx, pheno.name]), na.rm = TRUE)
    ## cat(m0, m1,m2, "\n")
    ret <- (m0 - m1) / (m0 - m2)
    if (!is.na(ret)) {
      if (ret > 1) { ret <- 1}
      if (ret < 0) { ret <- 0}
    }
    ret
  })
  head(ret)

  if (plot.it) {
    snapshot("plot.it", "plot.it.Rdata")
    ## draw linkage plot
    linkage.plot.pdf <- file.path(output, paste(prefix, "linkage_plot.pdf", sep = "."))
    pdf(file = linkage.plot.pdf, height = 8, width = 20)
    plot.manhattan(data.frame(Chrom = ret$chr, Position = ret$pos, Gene = ret$Gene, Pval = ret$additive), main = "additive")
    plot.manhattan(data.frame(Chrom = ret$chr, Position = ret$pos, Gene = ret$Gene, Pval = ret$recessive), main = "recessive")
    plot.manhattan(data.frame(Chrom = ret$chr, Position = ret$pos, Gene = ret$Gene, Pval = ret$dominant), main = "dominant")
    dev.off()
    mycat("Generated ", linkage.plot.pdf, "\n")

    dist.plot.pdf <- file.path(output, paste(prefix, "distribution_plot.pdf", sep = "."))
    ## require(gridExtra)
    ## pdf(file = dist.plot.pdf, height = 8, width = 8)
    tmp <- do.call(marrangeGrob, c(dist.plots, list(nrow=2, ncol=2)))
    ## print(ml)
    ##dev.off()
    ggsave(dist.plot.pdf, tmp, width = 6, height = 6)
    ## ggsave(dist.plot.pdf, ml, width = 8, height = 8)
    mycat("Generated ", dist.plot.pdf, "\n")
  }
  write.table(ret, file = fns$csv_file, quote = F, row.names = F, sep = ",")

  ## analysis <- list(test = test, detect = detect, tail = tail, prefix = prefix,
  ##                  bin = bin, results = results, bonferroni = bonferroni)
  wd <- getwd()
  save(list = ls(), file = fns$results_file)

  end.time <- Sys.time()
  diff.time <- difftime(end.time, start.time, units = "secs")
  msg <- (sprintf("meta_link() finished in %.3f seconds - [vcfFile=%s;pedFile=%s;pheno.name=%s;test=%s;detect=%s;tail=%s]",
                  diff.time,
                  vcfFile,
                  pedFile,
                  pheno.name,
                  test,
                  detect,
                  tail))
  print(msg)
  ## report("m", msg, fns$log_file)
  mycat(msg, "\n")
  mycat("Exit successfully\n")
  return(list(returncode = 0, message = "", result = ret))
}


## tmp <- pheno[, c("fid", "iid", "mother", "pheno", "gt")]
## tmp$fid <- factor(tmp$fid)
## tmp$mother <- factor(tmp$mother)
## g <- ggplot(data = tmp, aes(x = "gt", y = "pheno", col = "mother", pch = "fid")) + geom_point(position = "jitter")
## ggsave(filename = "tmp.pdf", plot = g, width = 7, height = 7)

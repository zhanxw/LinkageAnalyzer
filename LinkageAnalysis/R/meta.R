if (FALSE) {
  setwd("~/test.run/0530.meta")
  vcfFile <- "R0359-R0360.vcf"
  pedFile <- "R0359-R0360.ped"
  pheno.name <- "raw"
  output = "."
  test = "wG2"
  detect = "never"
  silent = T
  tail = "decreasing"
  prefix = ""
  plot.it = FALSE
  transform.pheno = NULL
}

#' Meta-analyze single variant in super pedigree
#' @exmaple
#' setwd("~zhanxw/test.run/0530.meta")
#' vcfFile <- "R0359-R0360.vcf"
#' pedFile <- "R0359-R0360.ped"
#' pheno.name <- "raw"
#' meta.single.link(vcfFile, pedFile, pheno.name)
meta.single.link <- function(vcfFile, ## a vector of list
                             pedFile, ## a vector of list
                             pheno.name, ## which phenotype to use
                             output = ".",
                             test = "wG2",
                             detect = "never",
                             silent = T,
                             tail = "decreasing",
                             prefix = "",
                             plot.it = FALSE,
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
crossCheckMotherGenotype <- function(pheno) {
  mom <- ddply(pheno, .(mother), function(x) {
    c(
        ref = sum(x$gt == 0, na.rm = TRUE),
        var = sum(x$gt > 0, na.rm = TRUE))
  })
  total.fix <- 0
  for (i in seq_len(nrow(mom))) {
    mom.name <- mom[i, "mother"]
    if (mom[i, "var"] > 0) {
      if (!mom.name %in% pheno$iid ) {
        next
      }
      if (pheno[ pheno$iid == mom.name, "gt"] != 1) {
        pheno[ pheno$iid == mom.name, "gt"] <- 1 ## set to het
        total.fix <- total.fix + 1
      }
    }
  }
  if (total.fix) {
    cat("Fixed ", total.fix, " mother genotypes. \n")
  }
  pheno
}

countForLethal <- function(pheno) {
  nHetMom <- nVarFromHetMom <- nUnknownMom <- nVarFromUnknownMom <- 0
  res <- list()
  for (i in seq_len(nrow(pheno))) {
    if (is.na(pheno[i, "gt"]) || pheno[i, "gt"] != 2) { ## we only care VAR
      next
    }
    mom.name <- pheno[i, "mother"]
    if (!mom.name %in% names(res) ) {
      res[[mom.name]] <- vector("character", 0)
    }
    res[[mom.name]] <- c(res[[mom.name]], pheno[i, "iid"])
  }
  for (i in seq_len(length(res))) {
    idx <- which(mom.name %in%  pheno$iid)
    if (length(idx) == 0) {
      # mom has unknown
      nUnknownMom <- nUnknownMom + 1
      nVarFromUnknownMom <- nVarFromUnknownMom + length(unique(res[[mom.name]]))
    } else {
      nHetMom <- nHetMom + 1
      nVarFromHetMom <- nVarFromHetMom + length(unique(res[[mom.name]]))
    }
  }
  list(nHetMom, nVarFromHetMom, nUnknownMom, nVarFromUnknownMom)
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
                                  plot.it = FALSE,
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
  mycat("Remove ", sum(idx), " samples from VCF as they are not in PED\n")
  vcf <- vcf.delete.sample.by.index(vcf, idx)

  idx <- match(vcf$sampleId, ped[,2])
  mycat("Remove ", nrow(ped) - length(idx), " samples from PED as they are not in VCF\n")
  ped <- ped[idx, ]
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


      msg <- sprintf("Log file [ %s ] created.", status.file.name)
      ## report("m", msg, log_file)
      mycat(msg)
      ##msg <- "Exit successfully but no outputs as dichotomization failed"
      msg <- "dichotomize failed"
      ## report("m", msg, log_file)
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

  # read data
  # get G3 mice
  g3.ped <- subset(ped, gen == 3)
  g3.name <- intersect(colnames(vcf$GT), g3.ped$iid)
  pheno <- subset(g3.ped, iid %in% g3.name )
  g2.name <- (pheno$mother)

  # for samples in PED but not in VCF, fill their genotype as NA
  tmp <- (setdiff(g3.name, colnames(vcf$GT)))
  if ( length(tmp)> 0) {
    mycat("WARNING: Genotypes of ", length(tmp), " G3 mice not in VCF:", tmp, "\n")
    vcf <- vcf.add.sample(vcf, tmp)
  }
  tmp <- (setdiff(g2.name, colnames(vcf$GT)))
  if ( length(tmp) > 0) {
    mycat("WARNING: Genotypes of ", length(tmp), " G2 mice not in VCF:", tmp, "\n")
    vcf <- vcf.add.sample(vcf, tmp)
  }

  # prepare genotype data
  geno <- vcf$GT[, g3.name]  # G3 genotypes
  ## geno.g2 <- vcf$GT[, g2.name]
  nVariant <- nrow(geno)

  numSex <- length(unique(pheno$sex))
  if (numSex == 1){
    null.model <- sprintf("%s ~ 1 + (1|fid) ", pheno.name)
  } else if (numSex >= 2) {
    null.model <- sprintf("%s ~ 1 + (1|fid) + sex", pheno.name)
  }

  if (test == "wG2") {
    if (all(is.na(pheno$mother))) {
      mycat("Does not have mother info at all (wG2 cannot work)!!\n")
      stop("Quit...")
    }
    null.model <- paste0(null.model, " + (1|mother)")
  }

  isBinary <- is.factor(pheno[,pheno.name])
  library(lme4)
  if (isBinary) {
    mycat("INFO: Phenotypes are treated as binary\n")
    null <- tryCatch(
        {
          glmer(as.formula(null.model), data = pheno, family = "binomial")
        },
        error = function(err) {
          print(str(err))
          print(err)
          msg <- ifelse(is.null(err[["message"]]), "UnknownError", err$message)
          ## msg <- paste("Exit failed", msg, sep = " ")
          return(list(returncode = 1, message = msg, error = err))
        })
    if (is.list(null) && null$returncode == 1) {
      mycat("Refit null model using less covariates\n")
      null.model <- str_replace(null.model, "\\+ sex", "")
      null <- glmer(as.formula(null.model), data = pheno, family = "binomial")
    }
  } else {
    mycat("INFO: Phenotypes are treated as continuous\n")
    null <- lmer(as.formula(null.model), data = pheno, REML = FALSE)
  }
  mycat("Finished fitting null model\n")

  snapshot("meta.link", "dbg.meta.fit.alt.Rdata")
  for (i in 1:nVariant) {
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
            error = function(err) {
              mycat("Fit [ ", alt.model, " binary=", isBinary, " ] failed ", str(err), "\n")
              ## print(err)
              msg <- ifelse(is.null(err[["message"]]), "UnknownError", err$message)
              return(list(returncode = 1, message = msg, error = err, isBinary = isBinary, alt.model = alt.model))
            })
        alt
      }
      alt <- run.mixed.effect.alt.model(isBinary, alt.model, pheno)

      if (! (is.list(alt) && alt$returncode == 1)) {
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

            ## check rank
            reduced.model.matrix <- model.matrix(as.formula(reduced.model), reduced.pheno)
            if ( qr(reduced.model.matrix)$rank < ncol(reduced.model.matrix)) {
              mycat("Insufficient sample size to fit models, set pvalue as one\n")
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
              mycat("Refit using glm() failed, set pvalue as one\n")
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

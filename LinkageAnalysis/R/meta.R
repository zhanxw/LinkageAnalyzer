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
    cat(..., file = log.file, append = TRUE)
  }

  wd <- getwd()
  mycat(paste("Version:", packageVersion("LinkageAnalysis")), "\n")
  mycat(paste("Date:", Sys.time()), "\n")
  mycat(paste("Host:", Sys.info()["nodename"]), "\n")
  mycat(paste("Call:", deparse(sys.status()$sys.calls[[1]])), "\n")
  mycat(paste("Directory:", wd), "\n")

  start.time <- Sys.time()
  mycat("Load VCF file: ", vcfFile, "\n")
  vcf <- get.vcf(vcfFile)
  vcf.summarize(vcf)

  mycat("Load PED file: ", pedFile, "\n")
  ped <- get.ped(pedFile)
  ped.summarize(ped)

  stopifnot(pheno.name %in% colnames(ped)[-(1:5)] )
  stopifnot(!all(is.na(ped[,pheno.name])))
  mycat("VCF/PED loaded\n")

  snapshot("meta.link", "dbg.meta.Rdata")

  mycat("Detect phenotype type\n")
  if (detect == "auto") {
    tmp <- dichotomize(ped[,pheno.name])
    if (tmp$succ) {
      ped[,pheno.name] <- tmp$new.value
    } else {
      mycat("ERROR: Dichotomize failed, ")

      status.file.name <- file.path(dirname(log.file),
                                    "R_jobs_complete_with_no_output.txt")
      cat(date(), file = status.file.name)
      cat("\t", file = status.file.name, append = TRUE)
      ## cat(msg, file = status.file.name, append = TRUE)
      ## cat("\n", file = status.file.name, append = TRUE)
      msg <- sprintf("Log file [ %s ] created.", status.file.name)
      ## report("m", msg, log_file)
      mycat(msg)
      msg <- "Exit successfully but no outputs as dichotomization failed"
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

  # read data
  # get G3 mice
  g3.ped <- subset(ped, gen == 3)
  g3.name <- intersect(colnames(vcf$GT), g3.ped$iid)
  pheno <- subset(g3.ped, iid %in% g3.name )
  g2.name <- (pheno$mother)

  # for samples in PED but not in VCF, fill their genotype as NA
  tmp <- (setdiff(g3.name, colnames(vcf$GT)))
  if ( length(tmp)> 0) {
    mycat("WARNING: Genotypes of ", length(tmp), " G3 mice not in VCF\n")
    vcf <- vcf.add.sample(vcf, tmp)
  }
  tmp <- (setdiff(g2.name, colnames(vcf$GT)))
  if ( length(tmp) > 0) {
    mycat("WARNING: Genotypes of ", length(tmp), " G2 mice not in VCF\n")
    vcf <- vcf.add.sample(vcf, tmp)
  }

  # prepare genotype data
  geno <- vcf$GT[, g3.name]  # G3 genotypes
  geno.g2 <- vcf$GT[, g2.name]
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
    null <- glmer(as.formula(null.model), data = pheno, family = "binomial")
  } else {
    mycat("INFO: Phenotypes are treated as continuous\n")
    null <- lmer(as.formula(null.model), data = pheno, REML = FALSE)
  }

  for (i in 1:nVariant) {
    if (i > 10 && is.debug.mode()) {
      cat("DEBUG....")
      next
    }
    mycat("Process ", i, " th variant: ", gene[i], "\n")
    if (length(unique(geno[i,])) == 1) {
      ## skip mono site
      next
    }
    for (type in c("additive", "recessive", "dominant")) {
      pheno$gt <- convert_gt(geno[i,], type)
      ## impute missing genotypes
      pheno$gt[is.na(pheno$gt)] <- mean(pheno$gt, na.rm = TRUE)

      alt.model <- paste0(null.model, " + gt")
      if (isBinary) {
        alt <- glmer(as.formula(alt.model), data = pheno, family = "binomial")
      } else {
        alt <- lmer(as.formula(alt.model), data = pheno, REML = FALSE)
      }

      pval <- anova(null, alt)$"Pr(>Chisq)"[2]
      if (is.na(pval)) {
        pval <- 1
      }
      direction <- fixef(alt)["gt"] > 0  # TRUE means protective effect, FALSE means harmful effect (desired)
      if (!is.na(direction)) {
        pval <- convert_tail(direction, pval, tail)
      }
      ret[i, type] <- pval
    }
    ret[i, "TDT"] <- NA
  }
  ret$REF <- rowSums(geno == 0, na.rm = TRUE)
  ret$HET <- rowSums(geno == 1, na.rm = TRUE)
  ret$VAR <- rowSums(geno == 2, na.rm = TRUE)
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

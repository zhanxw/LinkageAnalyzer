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

## param <- data.frame(func = d$func, vcf = vcfFile, ped = pedFile, pheno.name = pheno.name, ouput = d$output, test = d$test, detect = d$detect)
## write.csv(param, file = "param.csv", quote = FALSE, row.names = FALSE)

#' Perform single variant super pedigree analysis
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
      status.file.name <- paste(dirname(log_file),
                                "R_jobs_complete_with_no_output.txt",
                                sep = .Platform$file.sep)
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
      q('no')
    }
  } else {
    ## phenotype is quantitative, but let's double check it's QTL
    tmp <- ped[,pheno.name]
    tmp <- tmp[!is.na(tmp)] ## remove NA
    tmp <- tmp[! tmp %in% c(-9, 0, 1, 2) ]
    if (length(unique(tmp)) == 0) {
      mycat("WARNING: phenotype seems bo be binary, enable BINARY mode now\n")
      ##detect <- "auto"
      ## todo
    }

  }
  ## mycat("Phenotypes are treated as continuous\n")
  ## mycat("INFO: Phenotypes are treated as binary\n")


  fns <- filename(output, prefix)  # generate output file names
  if (!tail %in% c("increasing", "decreasing", "both")) {
    mycat("ERROR: Unrecognized option for tail [ ", tail, "]\n")
    return(list(returncode = -1))
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
    null <- glmer(as.formula(null.model), data = pheno, family = "binomial")
  } else {
    null <- lmer(as.formula(null.model), data = pheno, REML = FALSE)
  }


  for (i in 1:nVariant) {
    if (i > 10 ) {
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


## # this function tests the combinatory effect of two genes
## meta.double.link <- function(main_file, G2_file = "", output = ".", test = "woG2",
##                              detect = "never",
##                              silent = TRUE, tail = "decreasing", prefix = "",
##                              cutoff_single = 0.01, n_trial = 1e+05, plot.it = TRUE,
##                              transform.pheno = NULL) {
##   ## ##debug
##   ## save(list = ls(), file = "double_link.Rdata")
##   ## load("0512/nlrp3_inflammasome.3971/double_link.Rdata", verbose = TRUE)
##   ## setwd("~/test.run/0512/nlrp3_inflammasome.3971")
##   ## q('no')
##   ## if (FALSE) {
##   ##   load("/home/zhanxw/test.run/perm/Amber.1/td_rsfv_bga.20140427/double_link.Rdata", verbose = TRUE)
##   ##   source("/home/zhanxw/test.run/LinkageAnalysis/R/get_data.R")
##   ##   library(mclust)
##   ##   library(lme4)
##   ##   main <- main_file
##   ##   G2 <- G2_file
##   ##   log_file <- fns$log_file
##   ## }

##   start.time <- Sys.time()
##   # input
##   fns <- filename(output, prefix)  # generate output file names
##   input <- get_data(main_file, G2_file, fns$log_file, detect, transform.pheno)$data  # read data,G2 is null is G2 dam genotype data are not available

##   # check validity of parameters
##   if (!is.numeric(cutoff_single) || cutoff_single > 1) {
##     report("e", "Error in cutoff_single!", fns$log_file)
##   }
##   if (!tail %in% c("increasing", "decreasing", "both")) {
##     report("e", "Unrecognized option for tail!", fns$log_file)
##   }

##   # initialize significance matrix
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
##         if (exists("last.warning", envir = baseenv()) && !is.null(last.warning)){
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
##           pval <- anova_test(data, input$bin, test, silent, fns$log_file, tail, null.model = null.model)$pvalue
##           ## cat('pval = ', pval , "\n")
##           sig[[type]][i, j] <- pval
##           sig[[type]][j, i] <- sig[[type]][i, j]
##         }
##       }

##       # test for synthetic lethality
##       sig[["lethal"]][i, j] <- double_lethal(data, input, i, j, n_trial)
##       sig[["lethal"]][j, i] <- sig[["lethal"]][i, j]
##     } ## end loop j
##   } ## end loop i

##   # run single_link to determine which genes to mask
##   if (silent == FALSE) {
##     report("m", "Running single_link", fns$log_file)
##   }
##   single_link(main_file = main_file, G2_file = G2_file, detect = detect, output = tempdir(),
##               silent = TRUE, test = test, tail = tail, plot.it = FALSE)  # run single_link
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
##   msg <- (sprintf("double_link() finished in %.3f seconds - [main=%s;G2=%s;test=%s;detect=%s;tail=%s]",
##                   diff.time,
##                   main_file,
##                   G2_file,
##                   test,
##                   detect,
##                   tail))
##   print(msg)
##   report("m", msg, fns$log_file)
## }

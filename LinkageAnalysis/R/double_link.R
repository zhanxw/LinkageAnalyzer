if (FALSE) {
  ## library(LinkageAnalysis)
  fn <- list.files("~/test.run/LinkageAnalysis/R/", pattern = ".*R$")
  fn <- paste0("~/test.run/LinkageAnalysis/R/", fn)
  for (f in fn) {
    source(f)
  }

  setwd('~/test.run')
  main_file = "~/test.run/Rpackage.FACS_screen_B1b_cells.R0511.test1/R0511_main_norm_continuous.csv"
  G2_file = "~/test.run/Rpackage.FACS_screen_B1b_cells.R0511.test1/R0511_G2dam.csv"
  output = "tmp"
  test = "woG2"
  silent = FALSE
  prefix = "tmp.prefix"
  double_link(main_file, G2_file, output, test, silent, prefix = prefix)

}

# this function tests the combinatory effect of two genes
double_link <- function(main_file, G2_file = "", output = ".", test = "woG2",
                        detect = "never",
                        silent = TRUE, tail = "decreasing", prefix = "",
                        cutoff_single = 0.01, plot.it = TRUE,
                        transform.pheno = NULL) {
  log.file <- filename(output, prefix)$log_file
  ret <- tryCatch(
      {
        ret <- double.link.impl(main_file, G2_file, output, test, detect, silent, tail,
                                prefix, cutoff_single, plot.it, transform.pheno)
        if (ret$returncode == 0) {
          msg <- paste("Exit successfully", ret$message, sep = " ")
        } else {
          if (ret$returncode == 1 && ret$message == "dichototomize failed") {
            # this is a special error,
            # meaning we will treat it as normal exit but no output files
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
        print(str(err))
        print(err)
        msg <- ifelse(is.null(err[["message"]]), "UnknownError", err$message)
        msg <- paste("Exit failed", msg, sep = " ")
        report("m", msg, log.file)
        return(list(returncode = 1, message = msg))
      })

  return(ret)
}

double.link.impl <- function(main_file, G2_file = "", output = ".", test = "woG2",
                        detect = "never",
                        silent = TRUE, tail = "decreasing", prefix = "",
                        cutoff_single = 0.01, plot.it = TRUE,
                        transform.pheno = NULL) {
  ##debug
  ## wd <- getwd()
  ## save(list = ls(), file = "dbg.double_link.Rdata")
  ## load("0512/nlrp3_inflammasome.3971/double_link.Rdata", verbose = TRUE)
  ## setwd("~/test.run/0512/nlrp3_inflammasome.3971")
  ## q('no')
  ## if (FALSE) {
  ##   load("/home/zhanxw/test.run/perm/Amber.1/td_rsfv_bga.20140427/double_link.Rdata", verbose = TRUE)
  ##   source("/home/zhanxw/test.run/LinkageAnalysis/R/get_data.R")
  ##   library(mclust)
  ##   library(lme4)
  ##   main <- main_file
  ##   G2 <- G2_file
  ##   log_file <- fns$log_file
  ## }

  start.time <- Sys.time()
  # input
  fns <- filename(output, prefix)  # generate output file names

  # read data,G2 is null is G2 dam genotype data are not available
  snapshot("double.link.impl", "dbg.before.load.Rdata")
  tmp <- get_data(main_file, G2_file, fns$log_file, detect, transform.pheno)
  if (tmp$returncode) {
    return(tmp)
  }
  input <- tmp$data
  report("m", "Load data complete", fns$log_file)

  # check validity of parameters
  if (!is.numeric(cutoff_single) || cutoff_single > 1) {
    report("e", "Error in cutoff_single!", fns$log_file)
  }
  if (!tail %in% c("increasing", "decreasing", "both")) {
    report("e", "Unrecognized option for tail!", fns$log_file)
  }

  # initialize significance matrix
  snapshot("double.link.impl", "dbg.before.sig.Rdata")
  signif <- matrix(data = 1, nrow = input$n, ncol = input$n)  # significance matrix
  colnames(signif) <- input$genes$Gene
  rownames(signif) <- input$genes$Gene
  sig <- list(recessive = signif, additive = signif, dominant = signif,
              inhibitory = signif,
              lethal = signif)

  # statistical test
  null.model <- NULL
  null.model.ok <- NULL
  report("m", paste("Total ", input$n, " gene(s) to test"), fns$log_file)
  for (i in 1:(input$n - 1)) {
    ## if (i != 1) {
    ##   cat("skip i = ", i , " in debug\n")
    ##   next
    ## }
    if (silent == FALSE) {
      report("m", paste("---", input$genes$Gene[i], "---"), fns$log_file)
    }
    gt1 <- unlist(input$genotype[i, ])  # genotype of the first gene

    for (j in (i + 1):input$n) {
      ## if (i > 10 || j > 10) {
      ##   cat ("DBG skip\n")
      ##   next
      ## }
      print(sprintf("%s - %s x %s - (%d, %d, %d) - %.3f%%",
                    Sys.time(),
                    input$genes$Gene[i], input$genes$Gene[j],
                    i, j, input$n,
                    100. * ((j - i - 1) + (input$n - 1 + (input$n - i + 1)) * (i - 1) / 2) /
                    (input$n * (input$n - 1)/2) ))
      gt2 <- unlist(input$genotype[j, ])  # genotype of the second gene
      tmp <- convert_gt(gt1, "additive") - convert_gt(gt2, "additive")
      if (sum(!is.na(tmp)) <= 10) {
        next
      }  # too many NA values

      # prepare table of input and response variables
      data <- data.frame(pt = input$phenotype$phenotype, sex = input$phenotype$sex,
                         mother = input$phenotype$mother, gt1 = convert_gt(gt1, "additive"),
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

      # if the two genes locate within 30MB, skipped
      if (input$genes$chr[i] == input$genes$chr[j] &&
          abs(input$genes$pos[i] - input$genes$pos[j]) < 3e+07) {
        next
      }

      # fit null model
      ## source("/home/zhanxw/test.run/LinkageAnalysis/R/anova_test.R")
      if (is.null(null.model)) {
        cat("fit null model\n")
        assign("last.warning", NULL, envir = baseenv())
        null.model <- anova_test(data, input$bin, test, silent, fns$log_file, tail, fit.null = TRUE)$null.model
        if (exists("last.warning", envir = baseenv()) && !is.null(last.warning)){
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
          system.time(anova_test(data, input$bin, test, silent, fns$log_file, tail, null.model = null.model))
          system.time(anova_test(data, input$bin, test, silent, fns$log_file, tail))
        }
        if (!null.model.ok)  {
          pval  <- 1
          # cat("null.ok = ", null.model.ok, "\n")
          ## cat('pval = ', pval , "\n")
        } else {
          pval <- anova_test(data, input$bin, test, silent, fns$log_file, tail, null.model = null.model)$pvalue
          ## cat('pval = ', pval , "\n")
          sig[[type]][i, j] <- pval
          sig[[type]][j, i] <- sig[[type]][i, j]
        }
      }

      # test for synthetic lethality
      if (FALSE) {
        cat('dbg lethal\n')
        wd <- getwd()
        save(list = ls(), file = "dbg.mid.double_lethal.Rdata")
      }
      sig[["lethal"]][i, j] <- double_lethal(data, input, i, j)
      sig[["lethal"]][j, i] <- sig[["lethal"]][i, j]
    } ## end loop j
  } ## end loop i

  # run single_link to determine which genes to mask
  if (silent == FALSE) {
    report("m", "Running single_link", fns$log_file)
  }
  single_link(main_file = main_file, G2_file = G2_file, detect = detect, output = tempdir(),
              silent = TRUE, test = test, tail = tail, plot.it = FALSE)  # run single_link
  SignifOfSingle <- read.csv(filename(tempdir(), prefix = "")$csv_file)  # get significance
  SignifOfSingle$inhibitory <- 1  # dummy for inhibitory mode, include all
  unlink(as.vector(filename(tempdir(), prefix = "")))

  # calculate bonferroni correction cutoff
  cutoff_pair <- 0.05/(input$n * (input$n - 1)/2)  # the cutoff is approximate using this rough estimation
  report("m", paste("The Bonferroni cutoff is", pretty_num(cutoff_pair)), fns$log_file)

  if (plot.it) {
    # display results in heatmap (only genes that are not masked)
    if (silent == FALSE) {
      report("m", "Drawing heatmap", fns$log_file)
    }
    ## par(mfrow = c(1, 1))
    pdf(file = fns$pdf_file)

    for (type in c("recessive", "additive", "dominant", "inhibitory", "lethal")) {
      # draw a heatmap of results
      mask <- SignifOfSingle[, type] > cutoff_single
      heatmap_pval(sig, type, test, cutoff_pair, mask, fns)
    }

    dev.off()

    # make distribution plots (only genes that are not masked)
    if (silent == FALSE) {
      report("m", "Drawing distribution plots", fns$log_file)
    }
    ## if (TRUE) {
    ##   wd <- getwd()
    ##   cat("DBG in ", wd, "\n")
    ##   save(list = ls(), file = "tmp.Rdata")

    ##   load("tmp.Rdata" ,verbose = TRUE)
    ##   cutoff_pair <- 0.9
    ##   sig[["recessive"]][1,2] <- 1e-14
    ##   sig[["dominant"]][1,2] <- 1e-14
    ##   sig[["inhibitory"]][1,2] <- 1e-14
    ## }
    pdf(file = fns$distrib_file)
    par(mfrow = c(3, 2))

    for (i in 1:(input$n - 1)) {
      for (j in (i + 1):input$n) {
        # whether to draw the distribution plot, which types pass the cutoff
        flag <- ""
        type.min <- NULL
        p.min <- NULL
        for (type in c("recessive", "additive", "dominant", "inhibitory", "lethal")) {
          if (sig[[type]][i, j] < cutoff_pair &&
              SignifOfSingle[i, type] > cutoff_single &&
              SignifOfSingle[j, type] > cutoff_single) {
            flag <- paste(flag, " ", substr(type, 1, 2), " (", pretty_num(sig[[type]][i, j]), ")", sep = "")
          }
          if (is.null(type.min) || p.min > sig[[type]][i, j]) {
            type.min <- type
            p.min <- sig[[type]][i, j]
          }
        }
        if (flag == "") {
          next
        }
        ## cat ("DBG flag = ", flag, "type.min", type.min, "\n")

        gt1 <- unlist(input$genotype[i, ])  # genotype of the first gene
        gt2 <- unlist(input$genotype[j, ])  # genotype of the second gene
        data <- data.frame(pt = input$phenotype$phenotype, sex = input$phenotype$sex,
                           mother = input$phenotype$mother, gt1 = convert_gt(gt1, "additive"),
                           gt2 = convert_gt(gt2, "additive"))
        data <- data[(!is.na(data$gt1)) & (!is.na(data$gt2)), ]

        # distribution of genotype vs. phenotype
        ## cat("DBG: draw double_distrib for type ", type, "\n")

        double_distrib(flag, data, input, type.min, input$genes[c(i, j), ])
      }
    }

    ## if (FALSE) {
    ##   pdf(file = fns$distrib_file)
    ##   par(mfrow = c(3, 2))
    ##   double_distrib(flag, data, input, type, input$genes[c(i, j), ])  # distribution of genotype vs. phenotype
    ##   plot(1:5, 1:5, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
    ##   dev.off()
    ## }
    plot(1:5, 1:5, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
    dev.off()
    ## par(mfrow = c(1, 1))
  }

  # write p val matrix to file (full matrix)
  if (silent == FALSE) {
    report("m", "Outputting p value matrix", fns$log_file)
  }
  for (type in c("recessive", "additive", "dominant", "inhibitory", "lethal")) {
    write.csv(sig[[type]], file = sub(pattern = "full", replacement = paste(type),
                               fns$csv_file), quote = FALSE)
  }
  end.time <- Sys.time()
  diff.time <- difftime(end.time, start.time, units = "secs")
  msg <- (sprintf("double_link() finished in %.3f seconds - [main=%s;G2=%s;test=%s;detect=%s;tail=%s]",
                  diff.time,
                  main_file,
                  G2_file,
                  test,
                  detect,
                  tail))
  print(msg)
  report("m", msg, fns$log_file)

  return(list(returncode = 0, result = sig))
}

if (FALSE) {
  main_file <- "/home/zhanxw/test.run/perm/Amber.1/td_rsfv_bga.20140427/R0071_td_rsfv_bga_main_norm_continuous.csv"
  G2_file <- "/home/zhanxw/test.run/perm/Amber.1/td_rsfv_bga.20140427/R0071_G2dam.csv"
  main_file2 <- "/home/zhanxw/test.run/perm/Amber.2/td_rsfv_bga.20140427/R0071_td_rsfv_bga_main_norm_continuous.csv"
  G2_file2 <- "/home/zhanxw/test.run/perm/Amber.2/td_rsfv_bga.20140427/R0071_G2dam.csv"
  output = "."
  test = "wG2"
  detect = "never"
  silent = F
  tail = "decreasing"
  prefix = ""
  n_trail = 10 # 1e5
  plot.it = FALSE
  transform.pheno = NULL
  library(mclust)
  library(lme4)
  for (i in list.files("/home/zhanxw/test.run/LinkageAnalysis/R", pattern = ".*.R$")) {
    source(i)
  }
  ret <- single_link(main_file, G2_file, output, test, detect, silent, plot.it = TRUE)
  ret2 <- single_link(main_file2, G2_file2, output, test, detect, silent, plot.it = TRUE)
}


single_link <- function(main_file = "", G2_file = "",
                        output = ".", test = "wG2",
                        detect = "never", silent = T, tail = "decreasing",
                        prefix = "", plot.it = TRUE,
                        transform.pheno = NULL) {
  log.file <- filename(output, prefix)$log_file
  ret <- tryCatch(
      {
        ret <- single.link.impl(main_file, G2_file, output, test, detect, silent, tail,
                                prefix, plot.it, transform.pheno)
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

single.link.impl <- function(main_file = "", G2_file = "",
                             output = ".", test = "wG2",
                             detect = "never", silent = T, tail = "decreasing",
                             prefix = "", plot.it = TRUE,
                             transform.pheno = NULL) {
  start.time <- Sys.time()

  fns <- filename(output, prefix)  # generate output file names

  report("m", paste("Version:", packageVersion("LinkageAnalysis")), fns$log_file)
  report("m", paste("Date:", Sys.time()), fns$log_file)
  report("m", paste("Host:", Sys.info()["nodename"]) , fns$log_file)
  report("m", paste("Call:", deparse(sys.status()$sys.calls[[1]])), fns$log_file)

  if (!tail %in% c("increasing", "decreasing", "both")) {
    report("e", "Unrecognized option for tail!", fns$log_file)
  }

  # read data
  tmp <- get_data(main_file, G2_file, fns$log_file, detect, transform.pheno)
  if (tmp$returncode) {
    return(tmp)
  }
  raw_data <- tmp$data
  report("m", "Load data complete", fns$log_file)
  snapshot("single.load.data", "debug.single.load.data")

  genes <- raw_data$genes
  phenotype <- raw_data$phenotype
  genotype <- raw_data$genotype
  bin <- raw_data$bin
  G2 <- raw_data$G2  # G2 is null is G2 dam genotype data are not available

  # initialize
  pt <- phenotype$phenotype
  sex <- phenotype$sex
  mother <- phenotype$mother
  effective <- rep(TRUE, dim(genes)[1])  # whether stat test is carried out on each gene?

  sig_gene <- as.data.frame(matrix(data = 1, ncol = 4, nrow = dim(genes)[1]))  # significance of genes are stored in this vector
  colnames(sig_gene) <- c("additive", "recessive", "dominant", "TDT")

  if (plot.it) {
    # statistical testing and draw the distribution plot at the same time
    pdf(file = fns$distrib_file, width = 6, height = 9)
    par(mfrow = c(3, 2), mar = c(4, 4, 4, 4))
  }

  for (i in 1:dim(genes)[1]) {
    # get genotype in character string
    ## ## debug
    ## if (i != 2) {
    ##   next
    ## }
    ## if (genes$Gene[i] == "Kndc1") {
    ##   browser()
    ## } else {
    ##   next
    ## }
    if (silent == F) {
      msg <- paste("---", genes$Gene[i], "[", i, "of", dim(genes)[1], "]", "---")
      report("m", msg, fns$log_file)
    }
    gt <- unlist(genotype[i, ])
    if (sum(!is.na(convert_gt(gt, "additive"))) < 10) {
      if (silent == F) {
        report("m", "Skip because there are too few valid values", fns$log_file)
      }
      effective[i] <- FALSE
      next
    } else {
      if (silent == F) {
        report("m", paste("# remaining mice", sum(!is.na(convert_gt(gt, "additive")))),
               fns$log_file)
      }
    }

    # glmer anova test of full model and reduced model for each of the three types
    for (type in c("additive", "recessive", "dominant")) {
      data <- data.frame(pt = pt, sex = sex, gt = convert_gt(gt, type), mother = mother)
      data <- data[!is.na(data$gt), ]
      if (length(unique(data$gt)) > 1) {
        pval <- anova_test(data, bin, test, silent, fns$log_file, tail)$pvalue
        sig_gene[i, type] <- pval
        if (sig_gene[i,type]==0) {
          sig_gene[i,type]=1e-300
        }
      }
    }

    # TDT analysis
    data <- data.frame(pt = pt, sex = sex, gt = convert_gt(gt, "additive"), mother = mother)
    data <- data[!is.na(data$gt), ]
    sig_gene[i, "TDT"] <- TDT(data, G2, genes[i, ], bin, fns$log_file)

    # distribution plot
    if (length(unique(data$gt)) > 1 && plot.it) {
      distrib_plot(data, genes[i, ], bin)
    }
  }

  if (plot.it) {
    dev.off()  #close the distribution plot
    cat("Distribution plot saved at: ", fns$distrib_file, "\n")
    ## par(mfrow = c(1, 1))
  }

  # flag lethal genes
  genes$REF <- apply(genotype, 1, function(x) sum(x == "REF"))
  genes$HET <- apply(genotype, 1, function(x) sum(x == "HET"))
  genes$VAR <- apply(genotype, 1, function(x) sum(x == "VAR"))

  if (silent == F) {
    report("m", "Flag lethal genes\n", fns$log_file)
  }
  genes$lethal <- single_lethal(genotype, genes, phenotype, G2)

  alpha <- 0.05
  bonferroni <- alpha / sum(effective)
  if (plot.it) {
    # draw manhattan plot and diagnostic plot
    if (silent == F) {
      report("m", "Plotting Manhattan plots\n", fns$log_file)
    }
    pdf(file = fns$pdf_file, height = 8, width = 16)
    par(mfrow = c(1, 1))

    diagnostic(phenotype, bin, fns$log_file, raw_data$unconverted)  # draw diagnostic plot
    # manhattan plot
    for (type in c("additive", "recessive", "dominant", "TDT")) {
      print(head(genes))
      ## browser()
      mht(genes, sig_gene, type, fns$log_file, test, bin, effective,
          genotype)
    }
    dev.off()
  }
  # calculate penetrance and semidominance
  genetics <- get_genetics(genes, phenotype, genotype, bin, raw_data$unconverted)

  # save results
  if (silent == F) {
    msg <- sprintf("Save results in CSV format - %s", fns$log_file)
    report("m", msg, fns$log_file)
  }
  results <- cbind(genes[, c("Gene", "chr", "pos", "REF", "HET", "VAR", "lethal")],
                   sig_gene, genetics)
  write.table(results, file = fns$csv_file, quote = F, row.names = F, sep = ",")

  analysis <- list(test = test, detect = detect, tail = tail, prefix = prefix,
                   bin = bin, results = results, bonferroni = bonferroni)
  save(analysis, file = fns$results_file)

  end.time <- Sys.time()
  diff.time <- difftime(end.time, start.time, units = "secs")
  msg <- (sprintf("single_link() finished in %.3f seconds - [main=%s;G2=%s;test=%s;detect=%s;tail=%s]",
                  diff.time,
                  main_file,
                  G2_file,
                  test,
                  detect,
                  tail))
  print(msg)
  report("m", msg, fns$log_file)

  return(list(returncode = 0, result = results))
}

single_link <- function(main_file = "", G2_file = "", output = ".", test = "wG2",
                        detect = "never", silent = T, tail = "decreasing",
                        prefix = "", n_trial = 1e+05, plot.it = TRUE,
                        transform.pheno = NULL) {
    fns <- filename(output, prefix)  # generate output file names
    if (!tail %in% c("increasing", "decreasing", "both")) {
        report("e", "Unrecognized option for tail!", fns$log_file)
    }

    # read data
    raw_data <- get_data(main_file, G2_file, fns$log_file, detect, transform.pheno)
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
        if (silent == F) {
            report("m", paste("---", genes$Gene[i], "---"), fns$log_file)
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
                sig_gene[i, type] <- anova_test(data, bin, test, silent, fns$log_file, tail)
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
        ## par(mfrow = c(1, 1))        
    }
    
    # flag lethal genes
    genes$REF <- apply(genotype, 1, function(x) sum(x == "REF"))
    genes$HET <- apply(genotype, 1, function(x) sum(x == "HET"))
    genes$VAR <- apply(genotype, 1, function(x) sum(x == "VAR"))

    if (silent == F) {
        report("m", "Flag lethal genes\n", fns$log_file)
    }
    genes$lethal <- single_lethal(genotype, genes, phenotype, G2, n_trial)

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
            mht(genes, sig_gene, type, fns$log_file, test, bin, effective,
                              genotype)
        }
        dev.off()
    }
    # calculate penetrance and semidominance
    genetics <- get_genetics(genes, phenotype, genotype, bin, raw_data$unconverted)

    # save results
    if (silent == F) {
        report("m", "Save results in CSV format", fns$log_file)
    }
    results <- cbind(genes[, c("Gene", "chr", "pos", "REF", "HET", "VAR", "lethal")],
                     sig_gene, genetics)
    write.table(results, file = fns$csv_file, quote = F, row.names = F, sep = ",")

    analysis <- list(test = test, detect = detect, tail = tail, prefix = prefix,
                     n_trial = n_trial, bin = bin, results = results, bonferroni = bonferroni)
    save(analysis, file = fns$results_file)
}

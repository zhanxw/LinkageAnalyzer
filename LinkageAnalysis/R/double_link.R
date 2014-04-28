# this function tests the combinatory effect of two genes
double_link <- function(main_file, G2_file = "", output = ".", test = "woG2",
                        detect = "never",
                        silent = TRUE, tail = "decreasing", prefix = "",
                        cutoff_single = 0.01, n_trial = 1e+05, plot.it = TRUE,
                        transform.pheno = NULL) {
    # input
    fns <- filename(output, prefix)  # generate output file names
    input <- get_data(main_file, G2_file, fns$log_file, detect, transform.pheno)  # read data,G2 is null is G2 dam genotype data are not available

    # check validity of parameters
    if (!is.numeric(cutoff_single) || cutoff_single > 1) {
        report("e", "Error in cutoff_single!", fns$log_file)
    }
    if (!tail %in% c("increasing", "decreasing", "both")) {
        report("e", "Unrecognized option for tail!", fns$log_file)
    }

    # initialize significance matrix
    signif <- matrix(data = 1, nrow = input$n, ncol = input$n)  # significance matrix
    colnames(signif) <- input$genes$Gene
    rownames(signif) <- input$genes$Gene
    sig <- list(recessive = signif, additive = signif, dominant = signif,
                inhibitory = signif,
                lethal = signif)

    # statistical test
    report("m", paste("Total ", input$n, " gene(s) to test"), fns$log_file)    
    for (i in 1:(input$n - 1)) {
        if (silent == FALSE) {
            report("m", paste("---", input$genes$Gene[i], "---"), fns$log_file)
        }
        gt1 <- unlist(input$genotype[i, ])  # genotype of the first gene

        for (j in (i + 1):input$n) {
            gt2 <- unlist(input$genotype[j, ])  # genotype of the second gene
            tmp <- convert_gt(gt1, "additive") - convert_gt(gt2, "additive")
            if (sum(!is.na(tmp)) <= 10) {
                next
            }  # too many NA values

            # prepare table of input and response variables
            data <- data.frame(pt = input$phenotype$phenotype, sex = input$phenotype$sex,
                               mother = input$phenotype$mother, gt1 = convert_gt(gt1, "additive"),
                               gt2 = convert_gt(gt2, "additive"))
            data <- data[(!is.na(data$gt1)) & (!is.na(data$gt2)), ]
            if (length(unique(data$gt1)) == 1 || length(unique(data$gt2)) == 1) {
                next
            }

            # if the two genes locate within 30MB, skipped
            if (input$genes$chr[i] == input$genes$chr[j] &&
                abs(input$genes$pos[i] - input$genes$pos[j]) < 3e+07) {
                next
            }

            # test for combinatory effect
            for (type in c("recessive", "additive", "dominant", "inhibitory")) {
                if (type == "recessive") {
                    data$gt <- (data$gt1 + data$gt2 >= 4) * 1
                }  # gt is the interaction term
                if (type == "additive") {
                    data$gt <- data$gt1 * data$gt2
                }
                if (type == "dominant") {
                    data$gt <- (data$gt1 * data$gt2 >= 1) * 1
                }
                if (type == "inhibitory") {
                    data$gt <- data$gt1 * (data$gt2 == 0) + data$gt2 * (data$gt1 == 0)
                }
                if (length(unique(data$gt)) == 1) {
                    next
                }  # no difference in predictor values

                sig[[type]][i, j] <- anova_test(data, input$bin, test, silent, fns$log_file, tail)
                sig[[type]][j, i] <- sig[[type]][i, j]
            }

            # test for synthetic lethality
            sig[["lethal"]][i, j] <- double_lethal(data, input, i, j, n_trial)
            sig[["lethal"]][j, i] <- sig[["lethal"]][i, j]
        }
    }

    # run single_link to determine which genes to mask
    if (silent == FALSE) {
        report("m", "Running single_link", fns$log_file)
    }
    single_link(main_file = main_file, G2_file = G2_file, detect = detect, output = tempdir(),
                silent = TRUE, test = test, tail = tail)  # run single_link
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
        pdf(file = fns$distrib_file)
        par(mfrow = c(3, 2))

        for (i in 1:(input$n - 1)) {
            for (j in (i + 1):input$n) {
                # whether to draw the distribution plot, which types pass the cutoff
                flag <- ""
                for (type in c("recessive", "additive", "dominant", "inhibitory", "lethal")) {
                    if (sig[[type]][i, j] < cutoff_pair && SignifOfSingle[i, type] >
                        cutoff_single && SignifOfSingle[j, type] > cutoff_single) {
                        flag <- paste(flag, " ", substr(type, 1, 2), " (", pretty_num(sig[[type]][i, j]), ")", sep = "")
                    }
                }
                if (flag == "") {
                    next
                }

                gt1 <- unlist(input$genotype[i, ])  # genotype of the first gene
                gt2 <- unlist(input$genotype[j, ])  # genotype of the second gene
                data <- data.frame(pt = input$phenotype$phenotype, sex = input$phenotype$sex,
                                   mother = input$phenotype$mother, gt1 = convert_gt(gt1, "additive"),
                                   gt2 = convert_gt(gt2, "additive"))
                data <- data[(!is.na(data$gt1)) & (!is.na(data$gt2)), ]

                double_distrib(flag, data, input, type, input$genes[c(i, j), ])  # distribution of genotype vs. phenotype
            }
        }

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
}

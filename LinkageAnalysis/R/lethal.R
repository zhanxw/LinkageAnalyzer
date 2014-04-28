# G1 female is REF, G1 male is HET 1) Without G2 data G2 female (offspring) is
# REF/HET with 50% probability G2 male is HET (the same as G1 male)
# REF:HET:VAR=3:4:1 2) if G2 mother is REF, REF:HET:VAR=1:1:0 3) if G2 mother is
# HET, REF:HET:VAR=1:2:1

# this function runs random sampling to find the number of G3 mice with VAR
# genotype
single_sample <- function(mother_gt, n_G3, n_trial) {
    if (is.null(mother_gt) || is.na(mother_gt)) {
        mother_gt <- "unknown"
    }  # corresponding G2 data does not exist
    if (mother_gt == "REF") {
        return(0)
    }
    if (mother_gt == "HET") {
        prob <- 1/4
    }
    if (!mother_gt %in% c("REF", "HET")) {
        prob <- 1/8
    }  # unknown, FAILED, FALSE
    return(rbinom(n = n_trial, size = n_G3, prob = prob))
}

# this function tests whether each gene is a homozygous lethal with G2 data the
# H0 is a gene is not homozygous lethal
single_lethal <- function(genotype, genes, phenotype, G2, n_trial) {
    lethal <- rep(1, dim(genotype)[1])
    mothers <- as.vector(unique(phenotype$mother))

    for (i in 1:dim(genotype)[1]) {
        # skip genes with too many invalid values or on chrX
        gt <- unlist(genotype[i, ])
        gt <- convert_gt(gt, "recessive")
        if (sum(!is.na(gt)) < 3 || genes$chr[i] == "X") {
            next
        }

        trials <- rep(0, n_trial)  # the total number of G3 VAR mice in all trials

        # random sample for each mother and aggregate
        for (mother in mothers) {
            mother_gt <- G2[paste(genes$Gene[i], genes$Coord[i]), mother]
            # in case G2 genotype is unknown, if any child is VAR, the mother is definitely
            # HET
            if (any(gt[phenotype$mother == mother & !is.na(gt)] == 2)) {
                mother_gt <- "HET"
            }
            mother_G3_VAR <- single_sample(mother_gt, sum(phenotype$mother == mother &
                                                          !is.na(gt)), n_trial)
            trials <- trials + mother_G3_VAR
        }

        lethal[i] <- sum(trials <= sum(gt[!is.na(gt)] == 2))/n_trial
    }

    return(lethal)
}

# this function runs random sampling to find the number of G3 mice with desired
# genotype
double_sample <- function(mother_gt1, mother_gt2, n_G3, n_trial) {
    if (is.null(mother_gt1) || is.na(mother_gt1) ||
        mother_gt1 %in% c("FALSE", "FAILED", "ERROR")) {
        mother_gt1 <- "unknown"
    }
    if (is.null(mother_gt2) || is.na(mother_gt2) ||
        mother_gt2 %in% c("FALSE", "FAILED", "ERROR")) {
        mother_gt2 <- "unknown"
    }

    if (mother_gt1 == "REF" && mother_gt2 == "REF") {
        # cannot give the genotype we desired for sure
        ##return(rep(0, n_trial))
        prob <- 0
    }  
    if (mother_gt1 == "HET" && mother_gt2 == "HET") {
        prob <- 5/16
    }
    if ((mother_gt1 == "HET" && mother_gt2 == "REF") ||
        (mother_gt1 == "REF" && mother_gt2 == "HET")) {
        prob <- 1/8
    }
    if (mother_gt1 == "unknown" && mother_gt2 == "unknown") {
        prob <- (1/8 + 1/8 + 5/16)/4
    }
    if ((mother_gt1 == "unknown" && mother_gt2 == "REF") ||
        (mother_gt1 == "REF" && mother_gt2 == "unknown")) {
        prob <- (1/8)/2
    }
    if ((mother_gt1 == "unknown" && mother_gt2 == "HET") ||
        (mother_gt1 == "HET" && mother_gt2 == "unknown")) {
        prob <- (1/8 + 5/16)/2
    }

    ## return(rbinom(n = n_trial, size = n_G3, prob = prob))
    return(prob)
}

# this function tests for synthetic lethality of two genes (VAR,HET; HET,VAR; and
# VAR,VAR)
double_lethal <- function(data, input, i, j, n_trial) {
    mothers <- as.vector(unique(data$mother))
    if (any(input$genes[c(i, j), "chr"] == "X")) {
        return(1)
    }  # don't handle gene on chrX for the moment

    ## MonteCarlo <- matrix(data = 0, nrow = n_trial, ncol = length(mothers))
    ## colnames(MonteCarlo) <- mothers
    prob <- rep(0, length(mothers))
    names(prob) <- mothers
    numG3 <- rep(0, length(mothers))
    names(numG3) <- mothers
    
    for (mother in mothers) {
        mother_gt1 <- input$G2[paste(input$genes$Gene[i], input$genes$Coordination[i]),
                               mother]
        mother_gt2 <- input$G2[paste(input$genes$Gene[j], input$genes$Coordination[j]),
                               mother]
        n_G3 <- sum(data$mother == mother)

        # in case G2 genotype is unknown, if any child is VAR, the mother is definitely
        # HET
        if (any(data$gt1[data$mother == mother & !is.na(data$gt1)] == 2)) {
            mother_gt1 <- "HET"
        }
        if (any(data$gt2[data$mother == mother & !is.na(data$gt2)] == 2)) {
            mother_gt2 <- "HET"
        }
        ## MonteCarlo[, mother] <- double_sample(mother_gt1, mother_gt2, n_G3, n_trial)
        ## system.time(MonteCarlo[, mother] <- double_sample(mother_gt1, mother_gt2, n_G3, n_trial))
        prob[mother]  <- double_sample(mother_gt1, mother_gt2, n_G3, n_trial)
        numG3[mother] <- n_G3
    }

    obs <- sum((data$gt1 + data$gt2) >= 3)
    ## pval <- sum(apply(MonteCarlo, 1, sum) <= obs)/n_trial
    pval <- doubleSample(prob, numG3, obs, 1000000)
    return(pval)
}

# this function reads the main data and the optional G2 genotype data if G2
# genotype data, the returned list will have an additional element
get_data <- function(main = "", G2 = "", log_file, detect, transform.pheno = NULL) {
    data <- get_main(main, log_file, detect, transform.pheno)
    if (G2 != "") {
        data <- c(data, list(G2 = get_G2(G2, log_file)))
    }

    data$n <- dim(data$genes)[1]  # number of genes
    data$obs <- dim(data$genotype)[2]  # number of mice

    return(data)
}

get_G2 <- function(file = "", log_file) {
    # read G2 genotype data
    raw <- read.csv(file, header = T, stringsAsFactor = F)
    raw <- raw[raw[, 1] != "", ]
    raw <- raw[, raw[1, ] != "" & (!apply(raw, 2, function(x) {
        all(is.na(x))
    }))]

    # check col names
    if (all(colnames(raw)[1:3] == c("Gene", "Coordination", "Amplicon")) == F) {
        report("w", "Column names are not correct!", log_file)
    }

    # check duplicate genes or G2 mice
    colnames(raw)[-c(1:3)] <- toupper(colnames(raw)[-c(1:3)])
    mice <- colnames(raw)[-c(1:3)]
    if (length(mice) != length(unique(mice))) {
        report("e", "Duplicate mouse names!", log_file)
    }

    gene_pos <- paste(raw[, 1], raw[, 2])
    if (length(gene_pos) != length(unique(gene_pos))) {
        report("e", "Duplicate gene found!", log_file)
    }

    ## no genotype info
    if (ncol(raw) == 3) {
        rownames(raw) <- gene_pos
        return(raw)
    }

    # check level
    raw[, -c(1:3)] <- convert_upper(raw[, -c(1:3)])  # convert to upper case
    genotypes <- unlist(raw[, -c(1:3)])
    genotypes <- names(table(genotypes))
    if (any(c("REF", "HET") %in% genotypes) == F) {
        report("e", "None of REF/HET is found!", log_file)
    }
    if ("VAR" %in% genotypes) {
        report("w", "G2 mothers with VAR are found!", log_file)
    }

    tmp <- raw[, -c(1:3)]  # convert VAR to ERROR, this is a temporary solution
    tmp[tmp == "VAR"] <- "ERROR"
    raw[, -c(1:3)] <- tmp

    # return G2 genotypes
    rownames(raw) <- gene_pos
    return(raw)
}

get_main <- function(file = "", log_file, detect, transform.pheno=NULL) {
    # read raw data and check format
    raw <- read.csv(file, header = F, stringsAsFactor = F)
    raw[1:3, 1:2] <- "na"
    raw <- raw[!raw[, 1] == "", ]
    raw <- raw[, raw[1, ] != "" & (!apply(raw, 2, function(x) {
        all(is.na(x))
    }))]
    raw <- raw[!apply(raw[, -c(1:3)], 1, function(x) all(x == "FALSE")), ]  # delete all false ones
    raw <- raw[, c(rep(TRUE, 3), raw[3,-seq(3)] != "na")] ## kick out missing phenotypes
    # check names
    if (raw[1, 3] != "G2 Dam eartag" || !grepl(pattern = "gender", raw[2, 3], ignore.case = TRUE) ||
        raw[3, 3] != "Phenotype" || raw[4, 1] != "Gene" || raw[4, 2] != "Coordination" ||
        raw[4, 3] != "Amplicon") {
        report("w", "Names of the mice or genes are not correct!", log_file)
    }

    # check duplicate gene names or mouse names
    mice <- unlist(raw[4, -c(1:3)])
    if (length(mice) != length(unique(mice))) {
        report("e", "Duplicate mouse names!", log_file)
    }

    gene_pos <- paste(raw[-c(1:4), 1], raw[-c(1:4), 2])
    if (length(gene_pos) != length(unique(gene_pos))) {
        report("e", "Duplicate gene found!", log_file)
    }

    # check level
    # determine binary (TRUE/FALSE) depends on status (phenotype) is numeric or not
    status <- unlist(raw[3, -c(1:3)])
    if (suppressWarnings(any(is.na(as.numeric(status))))) {
        bin <- TRUE  # binary phenotype variable
        status <- names(table(status))
        if (length(status) == 1) {
            report("e", "Only affected mice or unaffected mice!", log_file)
        }
        if (!all(status %in% c("Affected", "Unaffected", "affected", "unaffected"))) {
            report("w", "Status other than affected and unaffected found!", log_file)
        }
    } else {
        bin <- FALSE
    }

    sex <- unlist(raw[2, -c(1:3)])
    sex <- names(table(sex))
    if (!all(sex %in% c("F", "M", "f", "m"))) {
        report("w", "Sex other than F/M found!", log_file)
    }

    raw[-c(1:4), -c(1:3)] <- convert_upper(raw[-c(1:4), -c(1:3)])  # convert to uppper case
    genotype <- names(table(unlist(raw[-c(1:4), -c(1:3)])))
    if (!all(c("REF", "HET", "VAR") %in% genotype)) {
        report("e", "At least one of REF/HET/VAR not found!", log_file)
    }

    # split table
    gene <- raw[-c(1:4), 1:2]
    colnames(gene) <- c("Gene", "Coordination")
    rownames(gene) <- raw[-c(1:4), 3] ## Amplicon as gene name
    gene <- as.data.frame(gene)

    pheno <- t(raw[1:3, -c(1:3)])  # mother, sex and phenotype are all capitalized
    colnames(pheno) <- c("mother", "sex", "phenotype")
    rownames(pheno) <- raw[4, -c(1:3)] ## ear tags
    pheno <- as.data.frame(pheno)
    pheno$mother <- factor(toupper(pheno$mother))

    pheno$sex <- toupper(pheno$sex)
    pheno$sex[pheno$sex == "F"] <- "0"  # unknown sex are converted to 0.5
    pheno$sex[pheno$sex == "M"] <- "1"
    pheno$sex[pheno$sex != "0" & pheno$sex != "1"] <- "0.5"
    pheno$sex <- as.numeric(pheno$sex)

    unconverted <- NULL  # used to store old continuous variable in case it is converted to binary

    if (bin == TRUE) {
        # read phenotype; if necessary, convert to binary variable
        pheno$phenotype <- factor(toupper(pheno$phenotype),
                                  levels = c("AFFECTED", "UNAFFECTED"),
                                  labels = c("AFFECTED", "UNAFFECTED"))  # binary phenotype
    } else {
        pheno$phenotype <- as.numeric(as.vector(pheno$phenotype))  # continous phenotype

        # transform phenotype
        if (!is.null(transform.pheno)) {
            if (transform.pheno == "log" || transform.pheno == "log10") {
                if (any(pheno$phenotype <= 0)) {
                    print(pheno$phenotype)
                    report("w", "At least one phenotype is non-positive and cannot be log transformed - now filled with minimum positive values", log_file)
                    min.pheno <- min (pheno$phenotype [ pheno$phenotype > 0], na.rm = TRUE)
                    pheno$phenotype [ pheno$phenotype <= 0] <- min.pheno
                }

                if (transform.pheno == "log") {
                    pheno$phenotype <- log(pheno$phenotype)
                } else if (transform.pheno == "log10") {
                    pheno$phenotype <- log10(pheno$phenotype)
                } else {
                    msg <- sprintf("Unrecognized transformation: [ %s ]", transform.pheno)
                    report("e", msg, log_file)
                }
            }
            ## print(pheno$phenotype)
        }

        # convert to binary variable
        if (detect != "never") {
            report("m", "Detecting clustering of mice", log_file)
            mclust1 <- Mclust(pheno$phenotype, G = 1, modelNames = "E")
            mclust2 <- Mclust(pheno$phenotype, G = 2, modelNames = "E")
            mclust3 <- Mclust(pheno$phenotype, G = 3, modelNames = "E")

            if (detect == "always" || (detect == "auto" && mclust2$bic > mclust1$bic)) {
                break_index <- sum(mclust2$classification == 1)  # find break point
                tmp <- pheno$phenotype
                tmp <- tmp[order(tmp)]
                break_point <- (tmp[break_index] + tmp[break_index + 1])/2

                report("m", "Convert to binary phenotype variable", log_file)
                report("m", paste("Cutoff:", break_point), log_file)
                bin <- TRUE
                unconverted$phenotype <- pheno$phenotype  # store in unconverted
                unconverted$break_point <- break_point
                pheno$phenotype <- factor(pheno$phenotype < break_point, levels = c(TRUE,
                                                                             FALSE), labels = c("AFFECTED", "UNAFFECTED"))
            }

            if (mclust3$bic > mclust1$bic && mclust3$bic > mclust2$bic) {
                report("m", "There may exist 3 clusters of phenotype scores!", log_file)
            }
        }
    }

    genotype <- raw[-c(1:4), -c(1:3)]
    colnames(genotype) <- raw[4, -c(1:3)]
    rownames(genotype) <- raw[-c(1:4), 3]
    genotype <- as.data.frame(genotype)

    # get chr/pos and sort genes table
    gene$chr <- sub(pattern = "_.*", replacement = "", gene$Coordination, perl = TRUE)  # split Coordinates
    gene$pos <- sub(pattern = ".*_", replacement = "", gene$Coordination, perl = TRUE)
    gene$pos <- as.numeric(gene$pos)
    temp <- suppressWarnings({
        new_order <- order(as.numeric(gene$chr), gene$pos)
    })
    gene <- gene[new_order, ]  # sort
    genotype <- genotype[new_order, ]  # sort

    # return the three tables
    return(list(genes = gene, phenotype = pheno, genotype = genotype, bin = bin,
                unconverted = unconverted))
}

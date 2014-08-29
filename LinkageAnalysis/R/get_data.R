#' this function reads the main data and the optional G2 genotype data if G2
#' genotype data, the returned list will have an additional element
#'
#' @param main_file a CSV file listing phenotyping and genotypes of G3 mice
#' @param G2 a CSV file listing G2 mothers' genotypes
#' @param log_file specify log file
#' @param detect "auto" means detect binary phenotype, or "NULL" means do nothing
#' @param transform.pheno "log" apply log transformation, or NULL means do nothing
#'
#' @return a list containing:
#'   returncode: 0 when success
#'   message: optionally message
#'   data
#'     genes: data.frame (sorted by chromsomal positions)
#'       Gene: gene names, character vector
#'       Coordination: character vector, e.g. "1_123"
#'       chr: character vector
#'       pos: numeric vector
#'     phenotype: a list
#'       mother: factor vector, mother names
#'       sex: 0 = female, 1 = male, 0.5 = unknown
#'       phenotype: numeric vector
#'     genotype: character matrix with rowname = gene, colname = sampleName, e.g. "REF"
#'     bin: logical
#'     unconverted: list (may be NULL)
#'       phenotype: original phenotype
#'       break_point: break point
#'     G2: character matrix. After fixing genotypes, it's similar to G2 file
#'       Gene: char vec
#'       Coordination: char vec
#'       Amplicon: char vec
#'     n:  # of genes
#'     obs:  # of mice
get_data <- function(main_file = "", G2 = "", log_file, detect, transform.pheno = NULL) {
  snapshot("get_data", "get_data.Rdata")
  data <- get_main(main_file, log_file, detect, transform.pheno)
  if (data$returncode) {
    return(data)
  } else {
    data <- data$data
  }

  if (G2 != "") {
    G2 <- get_G2(G2, log_file)
    if (G2$returncode) {
      return(G2)
    } else {
      G2 <- G2$data
    }
    data <- c(data, list(G2))
  }

  data$n <- dim(data$genes)[1]  # number of genes
  data$obs <- dim(data$genotype)[2]  # number of mice

  return(list(returncode = 0, message = "", data = data))
}

#' function to read G2 mice genotypes
#' @param file a CSV file listing G2 mothers' genotypes
#' @param log_file specify log file
#'
#' @return a list of
#' returncode: non-zero when criitical error
#' message: related to the function running status
#' data: data loaded
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
    return(list(returncode = 0, message = "", data = raw))
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
  return(list(returncode = 0, message = "", data = raw))
}

#' this function reads the main data
#'
#' @param main_file a CSV file listing phenotyping and genotypes of G3 mice
#' @param log_file specify log file
#' @param detect "auto" means detect binary phenotype, or "NULL" means do nothing
#' @param transform.pheno "log" apply log transformation, or NULL means do nothing
#'
#' @return NULL when criitical error
get_main <- function(main_file, log_file, detect, transform.pheno=NULL) {
  # read raw data and check format
  raw <- read.csv(main_file, header = F, stringsAsFactor = F)
  if (ncol(raw) <= 3) {
    return(list(
        returncode = 1,
        message = "No genotypes in the input file for further analysis",
        raw = raw))
  }
  raw[1:3, 1:2] <- "na"
  raw <- raw[!raw[, 1] == "", ]
  raw <- raw[, raw[1, ] != "" & (!apply(raw, 2, function(x) {
    all(is.na(x))
  }))]
  raw <- raw[!apply(raw[, -c(1:3), drop = FALSE], 1, function(x) all(x == "FALSE")), ]  # delete all false ones
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
    if (detect != "never" && detect != "NULL") {
      ## require(mclust)
      report("m", "Detecting clustering of mice", log_file)
      try.cluster <- tryCatch( {
        mclust1 <- Mclust(pheno$phenotype, G = 1, modelNames = "E")
        mclust2 <- Mclust(pheno$phenotype, G = 2, modelNames = "E")
        mclust3 <- Mclust(pheno$phenotype, G = 3, modelNames = "E")
        (list(returncode = 0))
      }, error = function(err) {
        print(str(err))
        print(err)
        msg <- ifelse(is.null(err[["message"]]), "UnknownError", err$message)
        msg <- paste("Cluster detection failed ", msg, " with ", nrow(pheno), " samples.", sep = " ")
        report("m", msg, log_file)
        return(list(returncode = 1, message = msg, pheno = pheno, gene = gene))
      })
      if (is.list(try.cluster) && try.cluster$returncode != 0) {
        return(try.cluster)
      }
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
        pheno$phenotype <- factor(pheno$phenotype < break_point,
                                  levels = c(TRUE, FALSE),
                                  labels = c("AFFECTED", "UNAFFECTED"))
        ## if very unbalanced cutoff has chosen, gracefully quit
        report("m", paste0("Number AFFECTED = ", sum(pheno$phenotype == "AFFECTED"),
                           " UNAFFECTED = ", length(pheno$phenotype)), log_file)
        ratio <- sum(pheno$phenotype == "AFFECTED") / length(pheno$phenotype)
        if (ratio < 0.2 || ratio > 0.8) {
          msg <- sprintf("Dichotomizing phentoypes failed (affected ratio = %f)", ratio)
          report("m", msg, log_file)

          # write status.file
          status.file.name <- file.path(dirname(log_file),
                                        "R_jobs_complete_with_no_output.txt")
          cat(date(), file = status.file.name)
          cat("\t", file = status.file.name, append = TRUE)
          cat(msg, file = status.file.name, append = TRUE)
          cat("\n", file = status.file.name, append = TRUE)

          msg <- sprintf("Log file [ %s ] created.", status.file.name)
          report("m", msg, log_file)
          return(list(returncode = 1, message = "dichotomize failed", data = raw))
        }
      }

      if (mclust3$bic > mclust1$bic && mclust3$bic > mclust2$bic) {
        report("m", "There may exist 3 clusters of phenotype scores!", log_file)
      }
    }
  }

  genotype <- raw[-c(1:4), -c(1:3), drop = FALSE]
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

  return(list(
      returncode = 0,
      message = "",
      data = list(genes = gene, phenotype = pheno, genotype = genotype, bin = bin,
          unconverted = unconverted)))
}

#' Dichotomize a numeric vector into AFFECTED/UNAFFECTED
#'
#' @param x numeric vector
#' @param log_file  specify log file
#'
#' @return list(), group -> suggested group number
dichotomize <- function(x, log_file = NULL) {
  snapshot("dichotomize", "dichotomize.Rdata")
  ret <- list(succ = FALSE, group = NA, break_point = NA, new.value = NULL, old.value = x)

  orig.len <- length(x)
  na.idx <- is.na(x)
  x <- x[!na.idx]

  report("m", "Detecting clustering of mice", log_file)
  ## require(mclust)
  mclust1 <- Mclust(x, G = 1, modelNames = "E")
  mclust2 <- Mclust(x, G = 2, modelNames = "E")
  mclust3 <- Mclust(x, G = 3, modelNames = "E")


  if (mclust2$bic > mclust1$bic) {
    break_index <- sum(mclust2$classification == 1)  # find break point
    tmp <- x
    tmp <- tmp[order(tmp)]
    break_point <- (tmp[break_index] + tmp[break_index + 1])/2

    #report("m", "Convert to binary outcomevariable", log_file)
    report("m", paste("Cutoff:", break_point), log_file)
    bin <- TRUE
    ret$break_point <- break_point
    ret$new.value <- factor(x < break_point,
                            levels = c(TRUE, FALSE),
                            labels = c("AFFECTED", "UNAFFECTED"))

    ## if very unbalanced cutoff has chosen, gracefully quit
    num.affected <- sum(ret$new.value == "AFFECTED")
    num.unaffected <- sum(ret$new.value == "UNAFFECTED")
    ratio <- num.affected / (num.affected + num.unaffected)
    report("m", paste0("Number AFFECTED = ", num.affected,
                       " UNAFFECTED = ", num.unaffected), log_file)

    if (ratio < 0.2 || ratio > 0.8) {
      msg <- sprintf("Dichotomizing phentoypes failed (affected ratio = %f)", ratio)
      report("m", msg, log_file)
      ret$group <- 1
      ret$succ  <- FALSE
    } else {
      ret$group <- 2
      ret$succ  <- TRUE
    }
  } else {
    msg <- sprintf("Dichotomizing phentoypes failed (one cluster is preferred to two clusters)\n")
    report("m", msg, log_file)
    msg <- sprintf("Data = %s", paste(ret$old.value, collapse = ","))
    report("m", msg, log_file)
    ret$group <- 1
    ret$succ  <- FALSE
  }

  if (mclust3$bic > mclust1$bic && mclust3$bic > mclust2$bic) {
    report("m", "There may exist 3 clusters of phenotype scores!", log_file)
    ret$group <- 3
    ret$succ  <- FALSE
  }
  ## fill in missing values
  nMissing <- length(ret$old.value) - length(ret$new.value)
  if ( nMissing > 0 && !is.na(ret$break_point)) {
    ret$new.value <- factor(ret$old.value < ret$break_point,
                            levels = c(TRUE, FALSE),
                            labels = c("AFFECTED", "UNAFFECTED"))
  }
  ret
}

#' read a PED file
#'
#' @param pedFile file name in PED format
#' @param pheno phenotype column header
#' @param detect "auto" means detect binary phenotype, or "NULL" means do nothing
#'
#' @return a data frame that looks like PED file content
get.ped <- function(pedFile, pheno = NULL, detect = NULL) {
  options(stringsAsFactors = FALSE)
  ## read peds
  ped <- read.table(pedFile, header = TRUE, comment.char = "!")
  colnames(ped)[1] <- "fid"
  if (any(duplicated(ped[,1:2]))) {
    stop("Duplicated entry")
  }

  ## make phenotype numeric
  for (i in 6:ncol(ped)) {
    ## if (!all(is.numeric(ped[ped[,i] != ".", i]))) {
    ##   warning("Some non-numerical phenotypes are not coded as .\n")
    ## }
    suppressWarnings(ped[,i] <- as.numeric(ped[,i]))
  }

  ## check generation
  markGeneration <- function(ped) {
    g1.names <- unique(ped[,1])
    ped$gen <- NA
    ped$gen[ped$iid %in% g1.names] <- 1
    g0.names <- c(ped$father[ped$iid %in% g1.names], ped$mother[ped$iid %in% g1.names])
    ped$gen[ped$iid %in% g0.names] <- 0
    ped$gen[ped$father %in% g1.names ] <- 2 ## g2 or g3
    g23.names <- ped$iid[ped$gen == 2]
    ped$gen[ped$mother %in% g23.names] <- 3
    ped$gen[ped$father %in% g23.names] <- 3
    ped$gen[ped$father == "." & ped$mother == "." & is.na(ped$gen)] <- -1 ## founders
    ## maybe add G4
    ped$gen
  }
  ped$gen <- markGeneration(ped)
  if (all(is.na(ped$gen))) {
    stop("Generation of ped file cannot be inferred!\n")
  }

  if (is.null(pheno)) {
    ped$pheno <- ped[,6]
  } else {
    ped$pheno <- ped[, pheno]
  }

  if (is.null(detect) || detect == "never" || detect == "NULL") {
    # do nothing
  } else if (detect == "auto"){
    tmp <- dichotomize(ped$pheno)
    if (tmp$succ) {
      ped$pheno <- tmp$new.value
    } else {
      if (tmp$group == 1) {
        cat (" cannot detect, quitting... \n")
        return(NA)
      } else {
        ## give warning
        cat (" more than two groups are possible!\n")
      }
    }
  }
  ped
}

#' Order by chromosomal positions
#'
#' @param chrom chacter vector listing chromosome names, e.g. "1", "2", ...
#' @param pos  integer vector
#'
#' @return indice which can be used to sort
order.chrom.pos <- function(chrom, pos) {
  ord <- order(chrom, pos)
  idx.num <- !is.na(suppressWarnings(as.numeric(chrom)) )
  ord[idx.num] <- order(as.numeric(chrom[idx.num])) - 100
  ord <- order(ord)
  ord
}

#' Read a VCF file
#'
#' @param vcfFile file name in VCF format
#' @param log_file specify a log file
#'
#' @return a list of parsed VCF file content
get.vcf <- function(vcfFile, log_file) {
  options(stringsAsFactors = FALSE)
  ## read vcfs
  vcf <- readLines(vcfFile)
  vcf <- vcf[!grepl("^##", vcf)]
  ## require(stringr)
  if (!grepl("^#CHROM\tPOS", vcf[1])) {
    stop("VCF input does not have valid header: #CHROM\tPOS...")
  }
  hdr <- str_split(str_replace(vcf[1], "#CHROM", "CHROM"), "\t")[[1]]

  vcf <- vcf[!grepl("^#", vcf)]
  parseVcf <- function(vcf, hdr) {
    ret <- str_split(vcf, "\t")
    ret <- do.call(rbind, ret)
    colnames(ret) <- hdr
    ## check duplication
    tmp <- duplicated(paste(ret[,1], ret[,2]))
    if (any(tmp)) {
      cat(vcf[tmp], "\n")
      warning("Found duplicated entries in VCF")
    }
    ret <- ret[!tmp,]

    ## sort by chromsomal (CHROM, POS) positions
    tmp <- order.chrom.pos(ret[,1], ret[,2])
    ret <- ret[tmp,]

    ## remove QUAL != "PASS"
    idx <- (ret[, "FILTER"] == "PASS")
    if (sum(!idx)) {
      loginfo("INFO: %d variant site do not PASS filter", sum(!idx))
    }
    ret <- ret[idx, ]

    ## process individual genotype fields
    site <- ret[, 1:8]
    indv <- ret[, -(1:9)]
    if (length(unique(ret[,9])) != 1) {
      logwarn("Format are not consistent")
      warnings("Format are not consistent")
    }
    fmt <- str_split(ret[1,9], ":")[[1]]

    ret.indv <- list()
    tmp <- indv
    for (i in 1:length(fmt)) {
      ret.indv[[fmt[i]]] <- sub(":.*", "", tmp)
      tmp <- sub('[^:]*:', "", tmp)

      loginfo("Process FORMAT tag %s", fmt[i])
      if (fmt[i] == "GT") {
        func <- function(x) {
          if (x == "0/0") {
            return(0)
          } else if (x == "0/1") {
            return (1)
          } else if (x == "1/1") {
            return (2)
          } else if (x == "." || x == "./.")  {
            return (NA)
          } else {
            msg <- gettextf("Unrecognized genotype: [ %s ]", x)
            warning(msg)
            return (NA)
          }
        }
        ret.indv[[fmt[i]]] <- apply(ret.indv[[fmt[i]]], c(1,2), func)
      } else if (fmt[i] == "REF" || fmt[i] == "VAR"){
        suppressWarnings(storage.mode(ret.indv[[fmt[i]]]) <- "integer")
      } else {
        warning(sprintf("Skip individual format tag: %s", fmt[i]))
        suppressWarnings(storage.mode(ret.indv[[fmt[i]]]) <- "integer")
      }
      ## ret.indv[[fmt[i]]] <- apply(ret.indv[[fmt[i]]], c(1,2), func)
    }
    ret <- c(as.data.frame(site), ret.indv)
    ret <- c(ret, list(sampleId = hdr[-(1:9)]))
    # length(ret)
    ret
  }
  vcf <- parseVcf(vcf, hdr)
  vcf
} ## get.vcf


#' Convert 0,1,2,NA to "REF", "HET", "ALT", "FAIL"
#'
#' @param x numeric vector
#'
#' @return character vector
geno.from.012 <- function(x) {
  ret <- x
  ret[x==0] <- "REF"
  ret[x==1] <- "HET"
  ret[x==2] <- "VAR"
  ret[is.na(x)] <- "FAIL"
  ret
}


#' Add samples to vcf data and set their genotypes as missing (GT, REF, VAR)
#'
#' @param vcf vcf data
#' @param sampleName character vector
#'
#' @return vcf data
vcf.add.sample <- function(vcf, sampleName) {
  ret <- vcf
  n <- length(sampleName)

  tmp <- vcf$GT[,rep(1, n), drop = FALSE]
  tmp[] <- NA
  colnames(tmp) <- sampleName
  ret$GT <- cbind(vcf$GT, tmp)

  ## skip first REF (column 3), and use the
  ref.idx <- which (names(vcf)  == "REF")[-1]
  tmp <- vcf[[ref.idx]][,rep(1, n), drop = FALSE]
  tmp[] <- NA
  colnames(tmp) <- sampleName
  ret[[ref.idx]] <- cbind(vcf[[ref.idx]], tmp)

  tmp <- vcf$VAR[,rep(1, n), drop = FALSE]
  tmp[] <- NA
  colnames(tmp) <- sampleName
  ret$VAR <- cbind(vcf$VAR, tmp)

  ret$sampleId <- c(vcf$sampleId, sampleName)
  ret
}

#' Delete samples to vcf data
#'
#' @param vcf vcf data
#' @param index integer vector, delete those samples
#'
#' @return vcf data
vcf.delete.sample.by.index <- function(vcf, index) {
  n <- length(vcf)
  for (i in 1:n) {
    if (is.matrix(vcf[[i]])) {
      if (ncol(vcf[[i]]) == length(index)) {
        vcf[[i]] <- vcf[[i]][,!index]
      }
    } else if (is.vector(vcf[[i]])) {
      next
    }
  }

  if (length(vcf[["sampleId"]]) == length(index)) {
    vcf[["sampleId"]] <- vcf[[i]][!index]
  }

  vcf
}

#' Summarize vcf data
#'
#' @param vcf vcf data
#'
vcf.summarize <- function(vcf) {
  nvar <- length(vcf$CHROM)
  nsample <- ncol(vcf$GT)
  loginfo("VCF contains %d variants and %d samples", nvar, nsample)
}

#' Summarize ped data
#'
#' @param ped ped data
#'
#' @return NULL
ped.summarize <- function(ped) {
  fam <- sort(unique(ped$fid))
  cat("PED contain ", length(fam), " families: ", fam, "\n")

  indv <- sort(unique(ped$iid))
  cat("PED contain ", length(indv), " individuals: ", head(indv), "... \n")

  father <- mother <- NULL ## bypass CRAN check
  fdr <- sort(unique(subset(ped, father == "." & mother == ".")$iid))
  cat("PED contain ", length(fdr), " founders: ", head(fdr), "... \n")

  gen <- sort(unique(ped$gen))
  cat("PED contains ", length(gen), " generations:\n")
  for (g in gen) {
    n <- sum(ped$gen == g, na.rm = TRUE)
    cat("  Generation ", g, " has ", n, " mice\n")
  }
  n <- sum(is.na(ped$gen))
  cat("  Generation unknown ", n, " mice\n")

  cat("PED gender distribution:\n")
  print(table(ped$sex, exclude = NULL))

  for (i in 6:ncol(ped)) {
    cat("PED phenotype summary: ", colnames(ped)[i], "\n")
    print(summary(ped[,i]))
  }
}

load.vcf.ped <- function(vcfFile, pedFile, pheno.name) {
  loginfo("Load VCF file: %s ", vcfFile)
  vcf <- get.vcf(vcfFile)
  vcf.summarize(vcf)

  if (length(vcf$CHROM) < 1) {
    msg <- "VCF is empty"
    logerror(msg)
    return(returncode = 1, message = msg)
  }

  loginfo("Load PED file: %s", pedFile)
  ped <- get.ped(pedFile)
  ped.summarize(ped)

  ## verify phenotype name
  stopifnot(pheno.name %in% colnames(ped)[-(1:5)] )
  stopifnot(!all(is.na(ped[,pheno.name])))

  ## clean up samples in VCF and PED
  idx <- ! vcf$sampleId %in% ped[,2]
  loginfo("Remove %d samples from VCF as they are not in PED.", sum(idx))  ## some sample may not be screened
  vcf <- vcf.delete.sample.by.index(vcf, idx)

  tmp <- setdiff(ped$iid, vcf$sampleId)
  loginfo("Add %d samples to VCF according to PED file.", length(tmp))
  vcf <- vcf.add.sample(vcf, tmp)

  ## rearrange PED according to VCF
  stopifnot(all(sort(vcf$sampleId) == sort(ped$iid)))
  idx <- match(vcf$sampleId, ped$iid)
  ped <- ped[idx,]
  stopifnot(all(vcf$sampleId == ped$iid))

  loginfo("%d/%d samples loaded as VCF/PED", length(vcf$sampleId), nrow(ped))
  return(list(vcf = vcf, ped = ped))
}

## prepare genotype and phenotype for additive/recessive/dominant model
prepare.model.data <- function(vcf, ped, pheno.name) {
  snapshot("prepare.model.data", "prepare.model.data.Rdata")
  ## remove suspicious samples (all REFs)
  idx <- apply(vcf$GT, 2, function(x) {all(x[!is.na(x)] == 0)})
  bad.sample <- vcf$sampleId[idx]
  if (length(bad.sample) > 0 ) {
    msg <- sprintf("Remove %d samples from VCF as their genotypes are REFs only: %s", sum(idx), paste(bad.sample, collapse = ","))
    g3.bad.sample <- subset(ped, ped$iid %in% bad.sample & ped$gen == 3)$iid
    if (length(g3.bad.sample) > 0) {
      msg <- c(msg, sprintf("Among them, %d are G3: %s", length(g3.bad.sample), paste(g3.bad.sample, collapse = ",")))
    } else {
      msg <- c(msg, sprintf("Among them, 0 are G3"))
    }
  }  else {
    msg <- sprintf("Remove 0 samples from VCF as their genotypes are REFs only.")
  }
  loginfo(msg)
  vcf <- vcf.delete.sample.by.index(vcf, idx)

  ## for superpedigree, when family A has genotyped for a variant but not family B,
  ## we need to manually impute REF for family B for such variant
  for (i in seq_len(nrow(vcf$GT))) {
    tmp <- data.frame(geno = vcf$GT[i,], iid = vcf$sampleId)
    tmp <- join(tmp, ped[,c("fid", "iid")], by = "iid")
    tmp2 <- ddply(tmp, .(fid), function(x){c(allMissing = (all(is.na(x$geno))))})
    idx <- tmp$fid %in% subset(tmp2, allMissing == TRUE)$fid
    vcf$GT[i, idx] <- 0  ## impute as REF
  }
  # reduce data by:
  #  1. remove mice with missing phenotype
  #  2. select G3 mice
  # get G3 mice
  pheno <- ped[!is.na(ped[,pheno.name]), ]
  pheno <- pheno[!is.na(pheno$gen), ]
  pheno <- pheno[pheno$gen == 3, ]
  pheno <- pheno[pheno$iid %in% colnames(vcf$GT), ]

  # prepare genotype data
  geno <- vcf$GT[, pheno$iid]  # G3 genotypes
  stopifnot(all(colnames(geno) == pheno$iid))

  # check missing rate
  missing <- rowMeans(is.na(geno))
  idx <- missing > 0.5
  msg <- sprintf("%d variants have >0.5 missing rate", sum(idx))
  logwarn(msg)
  missing <- colMeans(is.na(geno))
  idx <- missing > 0.5
  msg <- sprintf("%d samples have >0.5 missing rate", sum(idx))
  logwarn(msg)

  # report # of mice for models
  msg <- sprintf("%d mice will be analyzed in models", nrow(vcf$GT))
  loginfo(msg)

  list(pheno = pheno, geno = geno)
}

# if @param is "auto" then dichotomize phenotype,
# otherwise verify the phenotype is QTL
process.phenotype <- function(ped, pheno.name, detect) {
  if (detect == "auto") {
    tmp <- dichotomize(ped[,pheno.name])
    if (tmp$succ) {
      ped[,pheno.name] <- tmp$new.value
    } else {
      logerror("ERROR: Dichotomize failed")
      # write status.file
      log.file <- get("file", envir = getHandler('writeToFile'))
      status.file.name <- file.path(dirname(log.file),
                                    "R_jobs_complete_with_no_output.txt")
      cat(date(), file = status.file.name)
      cat("\t", file = status.file.name, append = TRUE)
      # write log
      msg <- sprintf("Log file [ %s ] created.", status.file.name)
      loginfo(msg)
      msg <- "dichotomize failed"
      logerror(msg)
      return(list(returncode = 1, message = msg))
    }
  } else {
    ## phenotype is quantitative, but let's double check it's QTL
    tmp <- ped[,pheno.name]
    tmp <- tmp[!is.na(tmp)] ## remove NA
    tmp <- tmp[! tmp %in% c(-9, 0, 1, 2) ]
    if (length(unique(tmp)) == 0) {
      logwarn("WARNING: phenotype seems bo be binary and will apply binary models\n")
      tmp <- ped[,pheno.name]
      idx <- tmp %in% c(-9, 0)
      ped[idx, pheno.name] <- NA
      ped[,pheno.name] <- factor(ped[,pheno.name], levels = c(1, 2), labels = c("AFFECTED", "UNAFFECTED"))
    }
  }
  return(list(ped = ped))
}

## #' Convert VCF and PED data to Tao's original format
## #'
## #' @param vcf vcf data
## #' @param ped ped data
## #'
## #' @return a list conforming to Tao's format
## #' @seealso get_data
## make.tao.format <- function(vcf, ped) {
##   ret <- list()
##   gene <- str_replace(str_extract(vcf$INFO, "GENE=[^;]*"), "GENE=", "")
##   chrom <- vcf$CHROM
##   pos <-  as.integer(vcf$POS)
##   coord <- paste(chrom, pos, sep = "_")
##   ret$genes <- data.frame(Gene = gene, Coordination = coord,
##                           chr = chrom, pos = pos)

##   ped.g3 <- subset(ped, gen == 3)
##   names.g3 <- ped.g3$iid
##   ret$phenotype <- list(mother = factor(ped.g3$mother),
##                         sex = ifelse(ped.g3$sex == 1, 1, ifelse(ped.g3$sex == 2, 0, 0.5)),
##                         phenotype = ped.g3$norm)

##   geno.g3 <- vcf$GT[, names.g3]
##   rownames(geno.g3) <- coord
##   colnames(geno.g3) <- names.g3
##   ret$genotype <- geno.from.012(geno.g3)

##   ret$bin <- NA

##   ret$unconverted <- list(phenotype = ret$phenotype$phenotype,
##                           break_point = NA)

##   ret$G2 <- cbind(gene, coord)
##   colnames(ret$G2) <- c("Gene", "Coordination")
##   ped.g2 <- subset(ped, gen == 2)
##   names.g2 <- ped.g2$iid
##   geno.g2 <- vcf$GT[, names.g2]
##   ret$G2 <- cbind(ret$G2, geno.from.012(geno.g2))

##   ret$n <- nrow(geno.g3)
##   ret$obs <- ncol(geno.g3)
##   ret
## }


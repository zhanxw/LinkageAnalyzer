# @param gt, when it's character, this function converts to numerical values
# representation according to type FAILED, FALSE, and other types will be
# converted to NA
# @param gt, when it's numerics, then conversion is directly peformed
convert_gt <- function(gt, type) {
  if (!type %in% c("additive", "recessive", "dominant") ) {
    cat("Unrecognized type: ", type, "\n")
    stop("convert_gt() failed")
  }
  if (is.character(gt)) {
    if (type == "additive") {
      gt2num <- c(0, 1, 2)
    } else if (type == "recessive") {
      gt2num <- c(0, 0, 1)
    } else if (type == "dominant") {
      gt2num <- c(0, 1, 1)
    }

    gt[gt == "REF"] <- gt2num[1]
    gt[gt == "HET"] <- gt2num[2]
    gt[gt == "VAR"] <- gt2num[3]
    gt[!gt %in% c("0", "1", "2")] <- NA
    gt <- as.numeric(gt)
  }

  ## conversion
  if (any(gt > 2 | gt < 0, na.rm = TRUE)) {
    cat("GT = ", gt, "\n")
    stop("some gt is not between [0, 2]")
  }
  if (type == "additive") {
    # do nothing
  } else if (type == "recessive") {
    gt <- ifelse(gt > 1.5, 1, 0) ##  (1.5, 2] => 1, [0, 1.5] => 0
  } else if (type == "dominant") {
    gt <- ifelse(gt > 0.5, 1, 0) ##  (0.5, 2] => 1, [0, 0.5] => 0
  }
  return(gt)
}

#' this function converts a two-tailed p value to a one-tailed p value
#' @param direction is positive
#' @param pval pvalue
#' @param tail "decreasing", "increasing" or "both"
convert_tail <- function(direction, pval, tail) {
  # expected direction
  if (tail == "decreasing") {
    if (direction == T) {
      pval <- 1 - pval/2
    } else {
      pval <- pval/2
    }
    return(pval)
  } else if (tail == "increasing") {
    # unexpected direction
    if (direction == F) {
      pval <- 1 - pval/2
    } else {
      pval <- pval/2
    }
    return(pval)
  } else if (tail == "both") {
    return(pval)
  } else {
    warning(paste0("Unrecognized tail type = ", tail))
    return(pval)
  }
}

# this function converts the character strings to upper case in one data frame
convert_upper <- function(x) {
  if (typeof(x) == "character")
    {
      x <- data.frame(data = x)
    }  # if x is a vector, convert to dataframe
  for (i in 1:dim(x)[2]) {
    x[, i] <- toupper(x[, i])
  }
  return(x)
}

# this function generates messages, warnings and errors (and stop the program)
report <- function(type, word, log_file = NULL, eol = TRUE) {
  # whether to print return sign
  if (eol == TRUE) {
    suffix <- "\r\n"
  } else {
    suffix <- ""
  }

  # message
  if (type == "m") {
    if (!is.null(log_file)) {
      cat(paste(word, suffix), file = log_file, append = TRUE)
    }
    loginfo(word)
    cat(paste(word, suffix))
  } else if (type == "w") {
    # warnings
    if (!is.null(log_file)) {
      cat(paste(word, suffix), file = log_file, append = TRUE)
    }
    logwarn(word)
    warning(word, call. = FALSE)
  } else if (type == "e") {
    # errors
    if (!is.null(log_file)){
      cat(paste(word, suffix), file = log_file, append = TRUE)
    }
    logerror(word)
    stop(word, call. = FALSE)
  }
}

# this function is used to generate the output file names in the output folder
filename <- function(output, prefix) {
  dir.create(output, showWarnings = FALSE)

  if (prefix != "") {
    prefix <- paste(prefix, "_", sep = "")
  }  # if a prefix is given

  log_file <- file.path(output, paste(prefix, "log.txt", sep = ""))
  unlink(log_file)

  linkage_file <- file.path(output, paste(prefix, "linkage_plot.pdf", sep = ""))
  csv_file <- file.path(output, paste(prefix, "full_results.csv", sep = ""))
  distrib_file <- file.path(output, paste(prefix, "distribution_plot.pdf", sep = ""))
  results_file <- file.path(output, paste("results.RData", sep = ""))

  return(list(log_file = log_file, linkage_file = linkage_file, csv_file = csv_file,
              distrib_file = distrib_file, results_file = results_file))
}

# this function shows small numbers (p values) in scientific mode
pretty_num <- function(x) {
  if (x == 0) {
    return(0)
  }
  if (x > 0.001) {
    return(round(x, digits = 3))
  }
  d <- floor(-log10(x))
  return(round(x * 10^d, digits = 3)/10^d)
}

#' re-source all .R files
#'
#' @param dir directory name
source.all.file <- function(dir = ".") {
  fn <- list.files(path = dir, pattern = ".R$")
  fn <- normalizePath(file.path(dir, fn))
  for (i in fn) {
    source(i)
  }
 if (TRUE) {
   library(stringr)
   library(plyr)
   library(ggplot2)
   library(gplots)
   library(lme4)
 }
}
if (FALSE) {
  source.all.file("~/test.run/LinkageAnalysis/R")
}

is.debug.mode <- function() {
  debug <- nchar(Sys.getenv(x="DEBUG_LINKAGE_ANALYSIS")) > 0
  return (debug)
}

enable.debug.mode <- function() {
  Sys.setenv(DEBUG_LINKAGE_ANALYSIS = 1)
}
disabe.debug.mode <- function() {
  Sys.unsetenv("DEBUG_LINKAGE_ANALYSIS")
}

# force will dump if debug/non-debug mode
snapshot <- function(call.func.name, fn, force = FALSE) {
  if (interactive()) {
    if (file.exists(fn)) {
      load(fn, envir = parent.frame(2), verbose = TRUE)
    } else {
      cat(fn, " does not exists, skipped.\n")
    }
    return(0)
  }
  if (is.debug.mode() || force) {
    if (!grepl("\\.Rdata", fn)) {
      fn <- paste0(fn, ".Rdata")
    }

    wd <- getwd()
    cat("DEBUG: function call = ", call.func.name, "\n")
    cat("DEBUG: current variables saved to: ", fn, "\n")
    env <- parent.frame()
    save(list = ls(envir=env), file = fn, envir = env)
  }
  return(0)
}

isIgnorableError <- function(x) {
  ## print("x = "); print(x); print(str(x))
  ignorable.errors <- c("dichotomize failed",
                        "Response is constant - cannot fit the model",
                        "Cannot model G2 effect as there is only one G2 mouse.",
                        "Sample size is smaller than free model parameters",
                        "Does not have mother info at all (wG2 cannot work)!!")
  if (!is.null(x$message) && x$message %in% ignorable.errors) {
    loginfo("Ignore error: %s", x$message)
    return (TRUE)
  }
  return(FALSE)
}

isSuccess <- function(x) {
  if (is.list(x) && !is.null(x$returncode) && x$returncode != 0) {
    return (FALSE)
  }
  return(TRUE)
}

changeSuffix <- function(fn, old.suffix, new.suffix) {
  pattern <- paste0(old.suffix, "$")
  sub(pattern = pattern, replacement = new.suffix, x = fn)
}

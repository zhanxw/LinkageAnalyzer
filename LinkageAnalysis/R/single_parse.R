#' @title Parse the analysis results of multiple runs of single_link
#'
#' @description
#' This function loads the RData file under each folder for one
#' \code{single_link} analysis event and make a summary of the analysis results
#' for multiple runs. WARNING: under each folder specified by \code{runs}, there
#' can only be the output files of one run of \code{single_link}!
#'
#' @param dir The main directory where there are at least one sub-folder which
#' contain(s) the analysis results of \code{single_link}
#' @param runs A vector of character strings. The sub-folder name(s) directly
#' under \code{dir} containing the analysis results of \code{single_link} that
#' should be pooled and parsed. WARNING: under each folder specified by
#' \code{runs}, there can only be the output files of one run of
#' \code{single_link}!
#' @param output The output folder to put the summarized output
#' @param prefix Default is "". An optional character string to be attached to
#' output file names. This can create personalized names for each job.
#' @param silent Print intermediate messages to stdout if set to FALSE.
#' @details The analysis output from multiple runs of \code{single_link} will be
#' parsed and summarized into a .csv file to be put under \code{output}. A .txt
#' log file will also be generated.
#' @return No return value.
#' @seealso \code{\link{single_link}}
#' @export
#' @examples
#' dir=system.file("extdata/single",package="LinkageAnalysis")
#' single_parse(dir=dir,runs=c("raw_continuous_woG2",
#'                             "raw_continuous_wG2",
#'                             "autoBinary_norm_continuous_woG2"),output=dir)
single_parse <- function(dir = ".", runs = c(), output = ".", silent = FALSE, prefix = "") {
  ######### preparations ##################
  filenames <- filename(output, prefix)
  tests <- c("lethal", "additive", "recessive", "dominant")
  if (length(runs) == 0 || any(!runs %in% list.files(dir))) {
    report("e", "At least one run is not found in the directory!", filenames$log_file)
  }
  options(scipen = 100)

  # get all analysis results
  all_results <- get_all_runs(dir, runs, filenames, silent)
  n <- dim(all_results[[runs[1]]]$results)[1]

  ######## extract the p values #######
  if (silent == FALSE) {
    msg <- sprintf("Write to .CSV output file - %s", filenames$log_file)
    report("m", msg, filenames$log_file)
  }

  # print headers
  write("The parameters for each run\n", file = filenames$csv_file)
  summary <- matrix(data = "", nrow = 6, ncol = length(runs) + 3)
  summary[1:6, 3] <- c("Run", "G2", "Direction", "Binary", "#Sampling", "Bonferroni")

  # running parameters
  for (run in runs) {
    summary[1:6, 3 + which(run == runs)] <- c(run, all_results[[run]]$test, all_results[[run]]$tail,
                               all_results[[run]]$bin, all_results[[run]]$n_trial, pretty_num(all_results[[run]]$bonferroni))
  }
  write.table(summary, file = filenames$csv_file, append = TRUE, quote = FALSE,
              row.names = FALSE, col.names = FALSE, sep = ",")
  write("", file = filenames$csv_file, append = TRUE)

  # print results for each test
  for (test in tests) {
    extract <- rep(0, n)  # find the significant genes to show
    for (run in runs) {
      bonferroni <- all_results[[run]]$bonferroni
      extract <- extract + 1 * (all_results[[run]]$results[, test] < bonferroni)
    }

    write(paste("Summary of the test of", test, "\n"), file = filenames$csv_file,
          append = TRUE)
    if (sum(extract) == 0) {
      next
    }

    # gene, chr, pos
    summary <- matrix(data = "", nrow = sum(extract > 0), ncol = length(runs) +
                      3)
    summary[, 1] <- all_results[[run]]$results[extract > 0, "Gene"]
    summary[, 2] <- paste("chr", all_results[[run]]$results[extract > 0, "chr"],
                          sep = "")
    summary[, 3] <- all_results[[run]]$results[extract > 0, "pos"]

    # extract p values
    for (run in runs) {

      tmp <- all_results[[run]]$results[extract > 0, test]
      for (i in 1:length(tmp)) {
        tmp[i] <- pretty_num(tmp[i])
      }
      summary[, which(run == runs) + 3] <- tmp
    }

    write.table(summary, file = filenames$csv_file, append = TRUE, quote = FALSE,
                row.names = FALSE, col.names = FALSE, sep = ",")
    write("", file = filenames$csv_file, append = TRUE)
  }

  write(paste("LinkageAnalysis version:", packageVersion("LinkageAnalysis"), "\n"),
        file = filenames$csv_file, append = TRUE)
  write(paste("Analysis conducted at:", date(), "\n"), file = filenames$csv_file,
        append = TRUE)

  options(scipen = 0)
}

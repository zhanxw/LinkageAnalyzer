##  ====================================================================================================================================================
##  |  LinkageAnalysis License                                                                                                                         |
##  |  ----------------------------------------------------------------------------------------------------------------------------------------------  |
##  |  a.   Copyright ©2014, The University of Texas Southwestern Medical Center.  All rights reserved; and                                            |
##  |  b.   This software and any related documentation constitutes published and/or unpublished works and may contain valuable trade secrets and      |
##  |       proprietary information belonging to The University of Texas Southwestern Medical Center (UT SOUTHWESTERN).  None of the foregoing         |
##  |       material may be copied, duplicated or disclosed without the express written permission of UT SOUTHWESTERN.  IN NO EVENT SHALL UT           |
##  |       SOUTHWESTERN BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING   |
##  |       OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF UT SOUTHWESTERN HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.         |
##  |       UT SOUTHWESTERN SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND        |
##  |       FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". UT          |
##  |       SOUTHWESTERN HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.                                   |
##  |  c.   This software contains copyrighted materials from R-package, ggplot2, gplots, gridExtra, lme4, logging, mclust, plyr and stringr.          |
##  |       Corresponding terms and conditions apply.                                                                                                  |
##  ====================================================================================================================================================

##  ====================================================================================================================================================
##  |  This file is part of LinkageAnalysis.													       |
##  |																		       |
##  |  LinkageAnalysis is free software: you can redistribute it and/or modify									       |
##  |  it under the terms of the GNU General Public License as published by									       |
##  |  the Free Software Foundation, either version 3 of the License, or									       |
##  |  (at your option) any later version.													       |
##  |																		       |
##  |  LinkageAnalysis is distributed in the hope that it will be useful,									       |
##  |  but WITHOUT ANY WARRANTY; without even the implied warranty of										       |
##  |  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the										       |
##  |  GNU General Public License for more details.												       |
##  |																		       |
##  |  You should have received a copy of the GNU General Public License									       |
##  |  along with LinkageAnalysis.  If not, see <http://www.gnu.org/licenses/>.									       |
##  ====================================================================================================================================================


#' @title Parse the analysis results of multiple runs of single.link
#'
#' @description
#' This function loads the RData file under each folder for one
#' \code{single.link} analysis event and make a summary of the analysis results
#' for multiple runs. WARNING: under each folder specified by \code{runs}, there
#' can only be the output files of one run of \code{single.link}!
#'
#' @param dir The main directory where there are at least one sub-folder which
#' contain(s) the analysis results of \code{single.link}
#' @param runs A vector of character strings. The sub-folder name(s) directly
#' under \code{dir} containing the analysis results of \code{single.link} that
#' should be pooled and parsed. WARNING: under each folder specified by
#' \code{runs}, there can only be the output files of one run of
#' \code{single.link}!
#' @param output The output folder to put the summarized output
#' @param prefix Default is "". An optional character string to be attached to
#' output file names. This can create personalized names for each job.
#' @param silent Print intermediate messages to stdout if set to FALSE.
#' @details The analysis output from multiple runs of \code{single.link} will be
#' parsed and summarized into a .csv file to be put under \code{output}. A .txt
#' log file will also be generated.
#' @return No return value.
#' @seealso \code{\link{single.link}}
#' @keywords internal
#' @examples
#' \dontrun{
#'  dir <- system.file("extdata/single/output",package="LinkageAnalysis")
#'  single_parse(dir=dir,runs=c("raw_continuous_woG2",
#'                              "raw_continuous_wG2",
#'                              "autoBinary_norm_continuous_woG2"),output=dir)
#' }
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

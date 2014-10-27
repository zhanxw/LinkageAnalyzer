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


get_run <- function(dir, run, filenames) {
  # load file
  results_file <- file.path(dir, run, "results.RData")
  suppressWarnings(rm("analysis"))
  tryCatch(load(results_file), error = function(e) report("e", "Record doesn't exist!",
    filenames$log_file))
  if (length(get("analysis")) != 8) {
    report("e", paste("Incomplete record for ", run, "!", sep = ""), filenames$log_file)
  }
  return(get("analysis"))
}

get_all_runs <- function(dir, runs, filenames, silent) {
  # get each run
  all_results <- c()
  for (run in runs) {
    if (silent == FALSE) {
      report("m", paste("Parsing:", run), filenames$log_file)
    }
    all_results[[run]] <- get_run(dir, run, filenames)
  }

  # check consistency
  same_file <- all_results[[run]]$results[, 1:3]  # the last run

  # check each run
  for (run in runs) {
    each_file <- all_results[[run]]$results[, 1:3]
    if (dim(each_file)[1] != dim(same_file)[1] || !all(same_file == each_file)) {
      report("e", "The analyses are not conducted on the same genes!", filenames$log_file)
    }
  }

  return(all_results)
}

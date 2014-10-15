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
.onAttach <- function(libname, pkgname){
  packageStartupMessage(sprintf("Package 'LinkageAnalysis' version %s",
                                packageVersion("LinkageAnalysis")))
  ## print license
  license <- getLicense()
  packageStartupMessage(license)

}

hasNewVersion <- function() {
  ## get new version
  newVersionLink = "http://zhanxw.com/LinkageAnalysis/version"
  tmpFile <- tempfile()
  ret <- tryCatch({
    download.file(newVersionLink, tmpFile, quiet = TRUE)
    readLines(tmpFile)
  }, error = function(e) {NULL})
  unlink(tmpFile)

  if (!is.null(ret) && length(ret) == 2) {
    version <- ret[1]
    if (utils::packageVersion("LinkageAnalysis") < version) {
      if (length(ret) > 1) {
        packageStartupMessage(ret[2])
      } else {
        packageStartupMessage("Found new version of LinkageAnalaysis: ", ret)
      }
      return(TRUE)
    }
  }
  return(FALSE)
}

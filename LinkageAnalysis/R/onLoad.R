.onAttach <- function(libname, pkgname){
  packageStartupMessage(sprintf("Package 'LinkageAnalysis' version %s",
                                packageVersion("LinkageAnalysis")))
  ## print license
  license <- getLicense()
  packageStartupMessage(license)

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

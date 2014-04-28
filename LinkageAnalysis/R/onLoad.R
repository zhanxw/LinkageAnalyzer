.onAttach <- function(libname, pkgname){
    ## newVersionLink = "http://zhanxw.com:8080/seqminer/version"
    ## conn <- url(newVersionLink)
    ## ret <- tryCatch(readLines(conn, n = 2), error = function(e) {NULL})
    ## close(conn)

    ## if (!is.null(ret) && length(ret) == 2) {
    ##     version <- ret[1]
    ##     if (utils::packageVersion("seqminer") < version) {
    ##         if (length(ret) > 1) {
    ##             packageStartupMessage(ret[2])
    ##         } else {
    ##             packageStartupMessage("Found new version of seqminer: ", ret)
    ##         }
    ##     }
    ## }
    packageStartupMessage(sprintf("Package 'LinkageAnalysis' version %s", packageVersion("LinkageAnalysis")))
}

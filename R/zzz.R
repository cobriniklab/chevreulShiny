

.onLoad <- function(libname, pkgname) {
    ns = topenv()
    ns$bw_files = system.file("extdata", "bw_files.db", package = "chevreulShiny")
}

.onAttach <- function(libname, pkgname) {
    suppressPackageStartupMessages(library(clustree))
}

.onAttach <- function(...){
  options("cheapr.cores" = getOption("cheapr.cores", 1))
}
.onUnload <- function(libname, pkgname){
  options(cheapr.cores = NULL)
}

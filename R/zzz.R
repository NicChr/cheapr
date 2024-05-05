.onAttach <- function(...){
  options("cheapr.cores" = getOption("cheapr.cores", 1),
          "cheapr.digits" = getOption("cheapr.digits", 2))
}
.onUnload <- function(libname, pkgname){
  options(cheapr.cores = NULL,
          cheapr.digits = NULL)
}

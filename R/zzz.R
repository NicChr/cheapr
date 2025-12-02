.onAttach <- function(...){
  options("cheapr.cores" = getOption("cheapr.cores", 1),
          "cheapr.digits" = getOption("cheapr.digits", 2))
  fastplyr_pkg <- find.package("fastplyr", quiet = TRUE)
  if (length(fastplyr_pkg) > 0){
    fastplyr_version <- utils::packageVersion("fastplyr")
    if (fastplyr_version < package_version("0.9.91")){
      packageStartupMessage(
        "fastplyr version >= 0.9.91 is needed with this version of cheapr (>= 1.5.0),
        please install it using `install.packages('fastplyr')`"
      )
    }
  }
}
.onUnload <- function(libname, pkgname){
  options(cheapr.cores = NULL,
          cheapr.digits = NULL)
}

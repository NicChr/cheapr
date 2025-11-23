.onAttach <- function(...){
  options("cheapr.cores" = getOption("cheapr.cores", 1),
          "cheapr.digits" = getOption("cheapr.digits", 2))
}
.onLoad <- function(...){
  fastplyr_pkg <- find.package("fastplyr", quiet = TRUE)
  if (length(fastplyr_pkg) > 0){
    fastplyr_version <- utils::packageVersion("fastplyr")
    if (fastplyr_version < package_version("0.9.91")){
      stop(
        "fastplyr version >= 0.9.91 is needed with this version of cheapr,
        please install it using `install.packages('fastplyr')`"
      )
    }
  }
}
.onUnload <- function(libname, pkgname){
  options(cheapr.cores = NULL,
          cheapr.digits = NULL)
}

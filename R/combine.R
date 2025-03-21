#' A cheapr version of `c()`
#'
#' cheapr's version of `c()`. It is quite a bit faster for atomic vectors
#' and combines data frame rows instead of cols.
#'
#' @param ... Objects to combine.
#'
#' @returns
#' Combined objects.
#'
#' @export
cheapr_c <- function(...){
  .Call(`_cheapr_cpp_c`, list(...))
}

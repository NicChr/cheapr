#' Memory-efficient alternative to `which()`
#'
#' @description
#' Exactly the same as `which()` but more memory efficient.
#'
#' @param x A [logical] vector.
#' @param invert If `TRUE`, indices of values that are not `TRUE` are returned
#' (including `NA`). If `FALSE` (the default), only `TRUE` indices are returned.
#'
#' @returns
#' An unnamed integer vector.
#'
#' @details
#' This implementation is similar in speed to `which()`
#' but usually more memory efficient.
#'
#' @examples
#' library(cheapr)
#' library(bench)
#' x <- sample(c(TRUE, FALSE), 1e05, TRUE)
#' x[sample.int(1e05, round(1e05/3))] <- NA
#'
#' mark(which_(TRUE), which(TRUE))
#' mark(which_(FALSE), which(FALSE))
#' mark(which_(logical()), which(logical()))
#' mark(which_(x), which(x), iterations = 20)
#' mark(base = which(is.na(match(x, TRUE))),
#'      collapse = collapse::whichv(x, TRUE, invert = TRUE),
#'      cheapr = which_(x, invert = TRUE),
#'      iterations = 20)
#'
#' @export
which_ <- function(x, invert = FALSE){
  .Call(`_cheapr_cpp_which_`, x, invert)
}

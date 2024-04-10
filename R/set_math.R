#' Math operations by reference
#'
#' @description
#' These functions transform your variable by reference, with no copies being made.
#'
#' @param x A numeric vector.
#' @param digits Number of digits to round to.
#' @param base Logarithm base.
#'
#' @details
#' These functions are particularly useful for situations
#' where you have made a copy and then
#' wish to perform further operations without creating more copies. \cr
#' `NA` and `NaN` values are ignored and never updated. \cr
#' These functions will \bold{not work} on \bold{any} classed objects, meaning they
#' only work on standard integer and numeric vectors and matrices.
#'
#' @returns
#' The exact same object with no copy made, just transformed.
#'
#' @examples
#' library(cheapr)
#' library(bench)
#'
#' x <- rnorm(2e05)
#' options(cheapr.cores = 2)
#' mark(
#'   base = exp(log(abs(x))),
#'   cheapr = set_exp(set_log(set_abs(x)))
#' )
#' options(cheapr.cores = 1)
#'
#' @rdname set_math
#' @export
set_abs <- cpp_set_abs
#' @rdname set_math
#' @export
set_floor <- cpp_set_floor
#' @rdname set_math
#' @export
set_ceiling <- cpp_set_ceiling
#' @rdname set_math
#' @export
set_trunc <- cpp_set_trunc
#' @rdname set_math
#' @export
set_exp <- cpp_set_exp
#' @rdname set_math
#' @export
set_sqrt <- cpp_set_sqrt
#' @rdname set_math
#' @export
set_change_sign <- cpp_set_change_sign
#' @rdname set_math
#' @export
set_round <- function(x, digits = 0){
  cpp_set_round(x, digits)
}
#' @rdname set_math
#' @export
set_log <- function(x, base = exp(1)){
  cpp_set_log(x, base)
}

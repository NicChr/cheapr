#' Math operations by reference - \bold{Experimental}
#'
#' @description
#' These functions transform your variable by reference, with no copies being made.
#'
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @param digits Number of digits to round to.
#' @param base Logarithm base.
#'
#' @details
#' These functions are particularly useful for situations
#' where you have made a copy and then
#' wish to perform further operations without creating more copies. \cr
#' `NA` and `NaN` values are ignored though in some instances `NaN` values may
#' be replaced with `NA`.
#' These functions will \bold{not work} on \bold{any} classed objects, meaning they
#' only work on standard integer and numeric vectors and matrices. \cr
#'
#' ### When a copy has to be made
#'
#' A copy is only made in certain instances, e.g. when passing an integer vector
#' to `set_log()`. A warning will always be thrown in this instance alerting the user
#' to assign the output to an object because `x` has not been updated by reference. \cr
#' To ensure consistent and expected outputs, always assign the output to the same object, \cr e.g.
#' `x <- set_log(x)` (\bold{do this}) ✅ \cr
#' `set_log(x)` (\bold{don't do this}) ❌ \cr
#' `x2 <- set_log(x)` (Don't do this either) ❌ \cr
#' \cr No copy is made here unless x is an integer vector.
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
  cpp_set_log(x, as.double(base))
}
#' @rdname set_math
#' @export
set_pow <- cpp_set_pow
#' @rdname set_math
#' @export
set_add <- cpp_set_add
#' @rdname set_math
#' @export
set_subtract <- cpp_set_subtract
#' @rdname set_math
#' @export
set_multiply <- cpp_set_multiply
#' @rdname set_math
#' @export
set_divide <- cpp_set_divide


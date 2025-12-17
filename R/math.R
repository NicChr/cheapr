#' Parallelised math operations
#'
#' @name math
#'
#' @param x `[numeric(n)]` vector.
#' @param y `[numeric(n)]` vector.
#' @param digits `[numeric(n)]` - Number of digits to round to.
#' @param base `[numeric(n)]` - Logarithm base.
#'
#' @returns
#' A transformed integer or double vector.
#'
#' @rdname math
#' @export
abs_ <- cpp_abs
#' @rdname math
#' @export
floor_ <- cpp_floor
#' @rdname math
#' @export
ceiling_ <- cpp_ceiling
#' @rdname math
#' @export
trunc_ <- cpp_trunc
#' @rdname math
#' @export
negate_ <- cpp_negate
#' @rdname math
#' @export
exp_ <- cpp_exp
#' @rdname math
#' @export
sqrt_ <- cpp_sqrt
#' @rdname math
#' @export
sign_ <- cpp_int_sign
#' @rdname math
#' @export
log_ <- function(x, base = exp(1)){
  cpp_log(x, base)
}
#' @rdname math
#' @export
log10_ <- function(x){
  cpp_log(x, base = 10)
}
#' @rdname math
#' @export
round_ <- function(x, digits = 0){
  cpp_round(x, digits)
}
#' @rdname math
#' @export
add_ <- cpp_add
#' @rdname math
#' @export
subtract_ <- cpp_subtract
#' @rdname math
#' @export
multiply_ <- cpp_multiply
#' @rdname math
#' @export
divide_ <- cpp_divide

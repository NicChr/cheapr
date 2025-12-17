#' Parallelised math operations
#'
#' @param x A numeric vector.

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
log_ <- function(x, base = exp(1)){
  cpp_log(x, base)
}
#' @rdname math
#' @export
log10_ <- function(x){
  cpp_log(x, base = 10)
}

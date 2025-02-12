#' A fast and integer-based `sign()`
#'
#' @param x Integer or double vector.
#'
#' @returns
#' An integer vector denoting the sign, -1 for negatives, 1 for positives and 0
#' for when `x == 0`.
#'
#' @export
int_sign <- cpp_int_sign

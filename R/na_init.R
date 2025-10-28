#' Fast `NA` initialisation
#'
#' @param x A vector.
#' @param n Vector length to initialise.
#'
#' @returns
#' Initialises `NA` values of the same type as `x`.
#'
#' @seealso [rep_len_]
#'
#' @export
na_init <- function(x, n = 0L) {
  .Call(`_cheapr_cpp_na_init`, x, n)
}

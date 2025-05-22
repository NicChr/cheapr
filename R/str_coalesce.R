#' Coalesce character vectors
#'
#' @description
#' `str_coalesce()` find the first non empty string `""`.
#' This is particularly useful for assigning and fixing the
#' names of R objects.
#'
#' In this implementation, the empty string `""` has priority over
#' `NA` which means `NA` is only returned when all
#' values are `NA`, e.g. `str_coalesce(NA, NA)`.
#'
#' @param ... Character vectors to coalesce.
#' @param .args An alternative to `...` for easier programming with lists.
#'
#' @returns
#' A coalesced character vector of length corresponding to the recycled
#' size of supplied character vectors. See `?recycle` for details.
#'
#' @details
#' `str_coalesce(x, y)` is equivalent to
#' `if_else(x != "" & !is.na(x), x, y)`.
#'
#' @examples
#' library(cheapr)
#'
#' # Normal examples
#' str_coalesce("", "hello")
#' str_coalesce("", NA, "goodbye")
#'
#' # '' always preferred
#' str_coalesce("", NA)
#' str_coalesce(NA, "")
#'
#' # Unless there are only NAs
#' str_coalesce(NA, NA)
#'
#' # `str_coalesce` is vectorised
#'
#' x <- val_insert(letters, "", n = 10)
#' y <- val_insert(LETTERS, "", n = 10)
#'
#' str_coalesce(x, y)
#'
#' # Using `.args` instead of `do.call` is much more efficient
#' library(bench)
#' x <- cheapr_rep_len(list(letters), 10^3)
#'
#' mark(do.call(str_coalesce, x),
#'      str_coalesce(.args = x),
#'      iterations = 50)
#'
#' @export
str_coalesce <- function(..., .args = NULL){
  .Call(`_cheapr_cpp_str_coalesce`, .Call(`_cheapr_cpp_list_args`, list(...), .args))
}

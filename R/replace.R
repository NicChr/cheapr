#' Fast vector replacement, an alternative to `[<-`
#'
#' @name replace
#'
#' @param x A vector.
#' @param where `[integer(n)]` - Where to assign replacement values.
#' @param with Replacement values.
#' @param in_place `[logical(1)]` - Should assignment be done in-place
#' (no copies)? Default is `FALSE`. Please note that assignment will occur
#' in-place where possible.
#'
#' @returns
#' A vector whose values are
#' replaced with `with` at locations specified by `where`.
#'
#' @examples
#' library(cheapr)
#'
#' x <- set_round(seq_(-2, 2, by = 0.5))
#'
#' x |>
#'   replace_(1, with = 100) # Assign value 100 at location 1
#'
#' # Base R casts to `x` and replacement to a common type
#' `[<-`(x, x== 0, "42")
#'
#' # `assign_at` only casts replacement to type of x
#' x |>
#'   replace_(x == 0, with = "42") # Assign value 42 where x == 0
#'
#' @rdname replace
#' @export
replace_ <- function(x, where, with, in_place = FALSE){
  .Call(`_cheapr_cpp_replace`, x, where, with, in_place)
}

base_assign_at <- function(x, where, with){
  x[where] <- with
  x
}

#' Very fast check that numeric vector consists only of whole numbers
#'
#' @param x `[numeric(n)]` - A numeric vector.
#' @param tol `[numeric(1)]` - Tolerance.
#' @param na.rm `[logical(1)]` - Should `NA` values be ignored?
#' Default is `TRUE`.
#'
#' @returns
#' `TRUE`, `FALSE`, or `NA` (see Details)
#'
#' @details
#' `is_whole_number()` will return `NA` when these 3 conditions are met:
#' - `na.rm` is `FALSE`
#' - `x` contains at least 1 `NA` value
#' - `x` contains only a mix of whole numbers and/or `NA` values. If any values
#' are not whole numbers then we can return `FALSE` even with the presence of
#' `NA` values.
#'
#' If `x` is not numeric then `is_whole_number()` always returns `FALSE`.
#'
#' @export
is_whole_number <- function(x, tol = sqrt(.Machine$double.eps), na.rm = TRUE){
  .Call(`_cheapr_cpp_is_whole_number`, x, tol, na.rm)
}

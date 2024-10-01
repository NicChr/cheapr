#' A sometimes cheaper but argument richer alternative to `.bincode()`
#'
#' @description
#' When `x` is an integer vector, `bin()` is cheaper than `.bincode()` as
#' no coercion to a double vector occurs. This alternative also has more
#' arguments that allow you to return the start values of the binned vector,
#' as well as including out-of-bounds intervals.
#'
#' @param x A numeric vector.
#' @param breaks A numeric vector of breaks.
#' @param left_closed Should intervals be left-closed (and right-open)?
#' Default is `TRUE`. If `FALSE` they are left-open (and right-closed).
#' @param include_endpoint Equivalent to `include.lowest` in  `?.bincode`.
#' @param include_oob Should out-of-bounds interval be included?
#' Default is `FALSE`. This is the equivalent of adding `Inf` as
#' the last value of the breaks, or `-Inf` as the
#' first value of the breaks if `left_closed = FALSE`.
#' @param codes Should an integer vector indicating which bin the values
#' fall into be returned? Default is `TRUE`. If `FALSE` the
#' start values of the respective bin intervals are returned, i.e the
#' corresponding breaks.
#'
#' @returns
#' Either an integer vector of codes indicating which bin the values fall
#' into, or the start of the intervals for which each value falls into.
#'
#' @export
bin <- function(x, breaks, left_closed = TRUE, include_endpoint = FALSE,
                include_oob = FALSE, codes = TRUE){
  cpp_bin(x, breaks, right = !left_closed, include_lowest = include_endpoint,
          include_oob = include_oob, codes = codes)
}

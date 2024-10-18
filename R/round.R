#' Simple wrappers to fix a niche problem
#'
#' @description
#' cheapr's `round` and `trunc` should operate exactly
#' like the base R equivalents. The constant `0` is added to get rid of
#' internal negative zeros in C that are carried to other C functions,
#' most significantly 'collapse' functions, producing erroneous results. \cr
#' See the examples for more details.
#'
#' An issue was raised here: github.com/SebKrantz/collapse/issues/452
#' but ultimately the issue was deemed too complex to solve.
#'
#'
#' @param x A numeric vector or complex vector for `round`.
#' @param digits Integer indicating the number of decimal places for `round`.
#' @param ... Arguments to be passed to methods.
#'
#' @returns
#' A numeric or complex vector.
#'
#' @examples
#' library(cheapr)
#' library(collapse)
#'
#' # There should be no visible difference to the user
#' # when using these functions.
#'
#' # The problem occurs when there exist negative and positive numbers
#' # that round to zero and then passed to functions that
#' # don't know how to different between C negative zeros and positive ones.
#'
#' x <- base::round(c(-0.05, 0.05))
#' x
#'
#' funique(x) # Res should trivially be 0
#' unique(x)
#'
#' group(x) # Res should be c(1, 1), not c(1, 2)
#'
#' data.table::rleid(x) # Res should be c(1, 1), not c(1, 2)
#'
#' # This issue should be eliminated with `cheapr::round()`
#'
#' y <- round(c(-0.05, 0.05))
#'
#' funique(y)
#' group(y)
#' data.table::rleid(y)
#'
#' @rdname round
#' @export
round <- function(x, digits = 0, ...){
  # add 0 to eliminate C representation of
  # negative zeros
  # which collapse::group() can't handle

  # To see the issue, run
  # collapse::funique(base::round(c(-0.05, 0.05)))
  base::round(x, digits, ...) + 0
}
#' @rdname round
#' @export
trunc <- function(x, ...){
  base::trunc(x, ...) + 0
}

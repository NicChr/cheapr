#' Cheaper version of `ifelse()`
#'
#' @name if_else
#'
#' @param condition [logical] A condition which will be used to
#' evaluate the if else operation.
#' @param true Value(s) to replace `TRUE` instances.
#' @param false Value(s) to replace `FALSE` instances.
#' @param na Catch-all value(s) to replace all other instances,
#' where `is.na(condition)`.
#'
#' @seealso [case] [val_match]
#'
#' @returns
#' A vector the same length as condition,
#' using a common type between `true`, `false` and `na`.
#'
#' @rdname if_else
#' @export
if_else_ <- function(condition, true, false, na = NULL){
  .Call(`_cheapr_cpp_if_else`, condition, true, false, na)
}
#' @rdname if_else
#' @export
cheapr_if_else <- if_else_


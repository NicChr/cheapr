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

# Internal function, do not use
# Assumes true, false, and na are similar types/classes
if_else2 <- function(condition, true, false, na){

  n <- length(condition)

  if (length(true) != 1 && length(true) != n){
    stop("`length(true)` must be 1 or `length(condition)`")
  }
  if (length(false) != 1 && length(false) != n){
    stop("`length(false)` must be 1 or `length(condition)`")
  }
  if (length(na) != 1 && length(na) != n){
    stop("`length(na)` must be 1 or `length(condition)`")
  }

  lgl_val_counts <- cpp_lgl_count(condition)
  n_true <- lgl_val_counts["true"]
  n_false <- lgl_val_counts["false"]
  n_na <- lgl_val_counts["na"]

  if (n_true == length(condition)){
    if (length(true) == 1){
      return(rep_len(true, length(condition)))
    } else {
      return(true)
    }
  }

  if (n_false == length(condition)){
    if (length(false) == 1){
      return(rep_len(false, length(condition)))
    } else {
      return(false)
    }
  }

  if (n_na == length(condition)){
    if (length(na) == 1){
      return(rep_len(na, length(condition)))
    } else {
      return(na)
    }
  }

  # The else part is most likely to be most prominent
  if (length(false) == length(condition)){
    out <- false
  } else {
    out <- rep_len(false, length(condition))
  }

  lgl_locs <- cpp_lgl_locs(condition, n_true = n_true, n_false = n_false,
                           include_true = TRUE, include_false = FALSE,
                           include_na = TRUE)
  true_locs <- lgl_locs[["true"]]
  na_locs <- lgl_locs[["na"]]

  if (length(true) == 1){
    out[true_locs] <- true
  } else {
    out[true_locs] <- true[true_locs]
  }
  if (length(na) == 1){
    out[na_locs] <- na
  } else {
    out[na_locs] <- na[na_locs]
  }
  out
}


#' Cheaper version of `ifelse()`
#'
#' @param condition [logical] A condition which will be used to
#' evaluate the if else operation.
#' @param true Value(s) to replace `TRUE` instances.
#' @param false Value(s) to replace `FALSE` instances.
#' @param default Catch-all value(s) to replace all other instances,
#' where `is.na(condition)`.
#'
#' @returns
#' A vector the same length as condition,
#' using a common type between `true`, `false` and `default`.
#'
#' @export
cheapr_if_else <- function(condition, true, false, default = false[NA_integer_]){

  if (!is.logical(condition)){
    stop("condition must be a logical vector")
  }
  if (length(true) != 1 && length(true) != length(condition)){
    stop("`length(true)` must be 1 or `length(condition)`")
  }
  if (length(false) != 1 && length(false) != length(condition)){
    stop("`length(false)` must be 1 or `length(condition)`")
  }
  if (length(default) != 1 && length(default) != length(condition)){
    stop("`length(default)` must be 1 or `length(condition)`")
  }

  if (is.factor(true) || is.factor(false) || is.factor(default)){
    template <- combine_factors(true[1L], false[1L], default[1L])[0L]
  } else {
    template <- c(true[1L], false[1L], default[1L])[0L]
  }

  true <- cast(true, template)
  false <- cast(false, template)
  default <- cast(default, template)

  if (is_base_atomic(true) && is_base_atomic(false) && is_base_atomic(default)){
    return(`mostattributes<-`(
      cpp_if_else(condition, true, false, default),
      attributes(template)
    ))
  }

  # Catch-all method

  lgl_val_counts <- cpp_lgl_count(condition)
  n_true <- lgl_val_counts["true"]
  n_false <- lgl_val_counts["false"]
  n_default <- lgl_val_counts["na"]

  if (n_true == length(condition)){
    if (length(true) == 1){
      return(rep(true, length(condition)))
    } else {
      return(true)
    }
  }

  if (n_false == length(condition)){
    if (length(false) == 1){
      return(rep(false, length(condition)))
    } else {
      return(false)
    }
  }

  if (n_default == length(condition)){
    if (length(default) == 1){
      return(rep(default, length(condition)))
    } else {
      return(default)
    }
  }

  # if (length(default) == 1 && is.na(default)){
  #   out <- rep(template, length.out = length(condition))
  # } else if (length(default) == length(condition)){
  #   out <- default
  # } else {
  #   out <- rep(default, length.out = length(condition))
  # }

  # if (length(default) == length(condition)){
  #   out <- default
  # } else {
  #   out <- rep(default, length.out = length(condition))
  # }

  # The else part is most likely to be most prominent
  if (length(false) == length(condition)){
    out <- false
  } else {
    out <- rep(false, length.out = length(condition))
  }

  lgl_locs <- cpp_lgl_locs(condition, n_true = n_true, n_false = n_false,
                           include_true = TRUE, include_false = FALSE,
                           include_na = TRUE)
  true_locs <- lgl_locs[["true"]]
  # false_locs <- lgl_locs[["false"]]
  default_locs <- lgl_locs[["na"]]

  if (length(true) == 1){
    out[true_locs] <- true
  } else {
    out[true_locs] <- true[true_locs]
  }
  # if (length(false) == 1){
  #   out[false_locs] <- false
  # } else {
  #   out[false_locs] <- false[false_locs]
  # }
  if (length(default) == 1){
    out[default_locs] <- default
  } else {
    out[default_locs] <- default[default_locs]
  }
  out
}

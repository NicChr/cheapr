#' Cheaper version of `ifelse()`
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
#' using a common type between `true`, `false` and `default`.
#'
#' @export
cheapr_if_else <- function(condition, true, false, na = false[NA_integer_]){

  if (!is.logical(condition)){
    stop("condition must be a logical vector")
  }
  if (length(true) != 1 && length(true) != length(condition)){
    stop("`length(true)` must be 1 or `length(condition)`")
  }
  if (length(false) != 1 && length(false) != length(condition)){
    stop("`length(false)` must be 1 or `length(condition)`")
  }
  if (length(na) != 1 && length(na) != length(condition)){
    stop("`length(na)` must be 1 or `length(condition)`")
  }

  if (is.factor(true) || is.factor(false) || is.factor(na)){
    template <- combine_factors(true[1L], false[1L], na[1L])
    template_lvls <- levels(template)
    true <- factor_(true, levels = template_lvls)
    false <- factor_(false, levels = template_lvls)
    na <- factor_(na, levels = template_lvls)
  } else {
    template <- c(true[1L], false[1L], na[1L])[0L]
    true <- cast(true, template)
    false <- cast(false, template)
    na <- cast(na, template)
  }


  if (is_base_atomic(true) && is_base_atomic(false) && is_base_atomic(na)){
    return(`mostattributes<-`(
      cpp_if_else(condition, true, false, na),
      attributes(template)
    ))
  }

  # Catch-all method

  lgl_val_counts <- cpp_lgl_count(condition)
  n_true <- lgl_val_counts["true"]
  n_false <- lgl_val_counts["false"]
  n_na <- lgl_val_counts["na"]

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

  if (n_na == length(condition)){
    if (length(na) == 1){
      return(rep(na, length(condition)))
    } else {
      return(na)
    }
  }

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

#' Count the number of occurrences of a value.
#'
#' @description
#' This is a programmer's version of `sum(x == value)` to count the number of
#' occurrences of a value without creating a potentially large logical vector.
#'
#' @param x A vector, list, data frame or matrix.
#' @param value A value with which to count the frequency of.
#' @param recursive Should the function be applied recursively to lists?
#' @param invert Should `which_val` find locations of
#' everything except specified value? Default is `FALSE`.
#'
#' @details
#' This is a generalisation of `num_na()` and as such the identity
#' `count_val(x, NA) == num_na(x)` will always hold.
#'
#' @returns
#' A count of the number of times `value` appears in `x`.
#'
#' @rdname scalars
#' @export
count_val <- function(x, value, recursive = TRUE){
  .Call(`_cheapr_cpp_count_val`, x, value, recursive)
}
#' @rdname scalars
#' @export
val_count <- count_val
#' @rdname scalars
#' @export
which_val <- function(x, value, invert = FALSE){
  .Call(`_cheapr_cpp_which_val`, x, value, invert)
}
val_replace <- function(x, value, replace){
  .Call(`_cheapr_cpp_val_replace`, x, value, replace, FALSE)
  # check_length(value, 1)
  # check_length(replace, 1)
  # which_replace <- which_val(x, value)
  # if (length(which_replace)){
  #   x[which_replace] <- replace
  # }
  # x
}
na_replace <- function(x, replace){
  val_replace(x, NA, replace)
}
val_rm <- function(x, value){
  n_vals <- count_val(x, value, recursive = TRUE)
  if (n_vals == unlisted_length(x)){
    sset(x, 0L)
  } else if (n_vals == 0){
    x
  } else {
    sset(x, cpp_which_val(x, value, invert = TRUE))
  }
}

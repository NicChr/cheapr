#' Efficient functions for counting, finding, replacing and removing scalars
#'
#' @description
#' These are primarily intended as very fast scalar-based functions
#' for developers.
#'
#' @param x A vector, list, data frame or matrix.
#' @param value A scalar value to count, find, replace or remove.
#' @param invert Should `which_val` find locations of
#' everything except specified value? Default is `FALSE`.
#' @param recursive Should values in a list be counted or replaced recursively?
#' Default is `TRUE` and very useful for data frames.
#' @param replace Replacement scalar value.
#'
#' @details
#' At the moment these functions only work for
#' integer, double and character vectors.
#' They are intended mainly for developers who wish to write cheaper code
#' and reduce expensive vector operations. For example
#' `val_count(x, 0)` will always be cheaper than `sum(x == 0)`.
#'
#' @returns
#' `val_count()` returns the number of times a scalar value appears in a vector
#' or list. \cr
#' `val_find()` returns the index locations of that scalar value. \cr
#' `val_replace()` replaces a specified scalar value with a replacement scalar
#' value. If no instances of said value are found then the input x is returned
#' as is. \cr
#' `na_replace()` is a convenience function
#' equivalent to `val_replace(x, NA, ...)`. \cr
#' `val_rm()` removes all instances of a specified scalar value.
#' If no instances are found, the original input x is returned as is.
#'
#' @rdname scalars
#' @export
val_count <- function(x, value, recursive = TRUE){
  .Call(`_cheapr_cpp_count_val`, x, value, recursive)
}
#' @rdname scalars
#' @export
count_val <- val_count
#' @rdname scalars
#' @export
val_find <- function(x, value, invert = FALSE){
  .Call(`_cheapr_cpp_which_val`, x, value, invert)
}
#' @rdname scalars
#' @export
which_val <- val_find
#' @rdname scalars
#' @export
val_replace <- function(x, value, replace, recursive = TRUE){
  .Call(`_cheapr_cpp_val_replace`, x, value, replace, recursive)
}
#' @rdname scalars
#' @export
na_replace <- function(x, replace, recursive = TRUE){
  val_replace(x, NA, replace, recursive = recursive)
}
#' @rdname scalars
#' @export
val_rm <- function(x, value){
  n_vals <- val_count(x, value, recursive = TRUE)
  if (n_vals == unlisted_length(x)){
    sset(x, 0L)
  } else if (n_vals == 0){
    x
  } else {
    sset(x, cpp_which_val(x, value, invert = TRUE))
  }
}

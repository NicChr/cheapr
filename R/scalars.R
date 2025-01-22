#' Efficient functions for counting, finding, replacing and removing scalars
#'
#' @description
#' These are primarily intended as very fast scalar-based functions
#' for developers.
#' They are particularly useful for working with `NA` values in a fast
#' and efficient manner.
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
#' The `val_` functions allow you to very efficiently work with
#' scalars, i.e length 1 vectors. Many common common operations like
#' counting the occurrence of `NA` or zeros, e.g. `sum(x == 0)` or
#' `sum(is.na(x))` can be replaced more efficiently with
#' `val_count(x, 0)` and `na_count(x)` respectively.
#'
#' At the moment these functions only work for
#' integer, double and character vectors with the exception of the `NA`
#' functions.
#' They are intended mainly for developers who wish to write cheaper code
#' and reduce expensive vector operations.
#'
#' * `val_count()` - Counts occurrences of a value
#' * `val_find()` Finds locations (indices) of a value
#' * `val_replace()` - Replaces value with another value
#' * `val_rm()` - Removes occurrences of value from an object
#'
#' There are `NA` equivalent convenience functions.
#'
#' * `na_count()` == `val_count(x, NA)`
#' * `na_find()` == `val_find(x, NA)`
#' * `na_replace()` == `val_replace(x, NA)`
#' * `na_rm()` == `val_rm(x, NA)`
#'
#' `val_count()` and `val_replace()` can work recursively. For example,
#' when applied to a data frame, `na_replace` will replace `NA` values across
#' the entire data frame with the specified replacement value.
#'
#' In 'cheapr' function-naming conventions have not been consistent but
#' going forward
#' all scalar functions (including the `NA` convenience functions) will be
#' prefixed with 'val_' and 'na_' respectively.
#' Functions named with the older naming scheme like `which_na` may be
#' removed at some point in the future.
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
val_rm <- cpp_val_remove
#' @rdname scalars
#' @export
na_count <- num_na
#' @rdname scalars
#' @export
na_find <- function(x, invert = FALSE){
  if (invert){
    which_not_na(x)
  } else {
    which_na(x)
  }
}
#' @rdname scalars
#' @export
na_rm <- function(x){

  ## This part is for legacy reasons
  ## Some users might be using `na_rm` to remove empty rows of a data frame
  if (is.list(x)){
    n_na <- na_count(x, recursive = TRUE)
    if (n_na == unlisted_length(x)){
      sset(x, 0L)
    } else if (n_na == 0){
      x
    } else {
      na_locs <- na_find(x, invert = TRUE)
      sset(x, na_locs)
    }
  } else {
    .Call(`_cheapr_cpp_val_remove`, x, NA)
  }

}

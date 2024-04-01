#' Efficient functions for dealing with missing values.
#'
#' @description
#' `is_na()` is a parallelised alternative to `is.na()`. \cr
#' `num_na(x)` is a faster and more efficient `sum(is.na(x))`. \cr
#' `which_na(x)` is a more efficient `which(is.na(x))` \cr
#' `which_not_na(x)` is a more efficient `which(!is.na(x))` \cr
#' `row_na_counts(x)` is a more efficient `rowSums(is.na(x))` \cr
#' `row_all_na()` returns a logical vector indicating which rows are empty
#' and have only `NA` values. \cr
#' `row_any_na()` returns a logical vector indicating which rows have at least
#' 1 `NA` value. \cr
#' The `col_` variants are the same, but operate by-column.
#'
#' @details
#' These functions are designed primarily for programmers, to increase the speed
#' and memory-efficiency of `NA` handling. \cr
#' Most of these functions can be parallelised through `options(cheapr.cores)`. \cr
#'
#' ### Common use-cases
#' To replicate `complete.cases(x)`, use `!row_any_na(x)`. \cr
#' To find rows with any empty values,
#' use `which_(row_any_na(df))`. \cr
#' To find empty rows use `which_(row_all_na(df))` or `which_na(df)`.
#' To drop empty rows use `na_rm(df)` or `sset(df, which_(row_all_na(df), TRUE))`.
#'
#' ### `is_na`
#' `is_na` Is an S3 generic function. It will internally fall back on
#' using `is.na` if it can't find a suitable method.
#' Alternatively you can write your own `is_na` method.
#' For example there is a method for `vctrs_rcrd`
#' objects that simply converts it to a data frame and then calls `row_all_na()`.
#' There is also a `POSIXlt` method for `is_na` that is much faster than `is.na`.
#'
#' ### Lists
#' When `x` is a list, `num_na`, `any_na` and `all_na` will recursively search
#' the list for `NA` values. If `recursive = F` then `is_na()` is used to
#' find `NA` values. \cr
#' `is_na` differs to `is.na` in 2 ways:
#' * List elements are counted as `NA` if either that value is `NA`, or
#' if it's a list, then all values of that list are `NA`.
#' * When called on a data frame, it returns `TRUE` for empty rows that contain
#' only `NA` values.
#'
#' @param x A vector, list, data frame or matrix.
#' @param recursive Should the function be applied recursively to lists?
#' The default is `TRUE`. Setting this to `TRUE` is actually much cheaper because
#' when `FALSE`, the other `NA` functions rely on calling `is_na()`,
#' therefore allocating a vector. This is so that alternative objects with
#' `is.na` methods can be supported.
#'
#'
#' @returns
#' Number or location of `NA` values.
#'
#' @examples
#' library(cheapr)
#' library(bench)
#'
#' x <- 1:10
#' x[c(1, 5, 10)] <- NA
#' num_na(x)
#' which_na(x)
#' which_not_na(x)
#'
#' row_nas <- row_na_counts(airquality)
#' col_nas <- col_na_counts(airquality)
#' names(row_nas) <- rownames(airquality)
#' names(col_nas) <- colnames(airquality)
#' row_nas
#' col_nas
#'
#' df <- airquality[, 1:2]
#'
#' # Number of NAs in data
#' num_na(df)
#' # Which rows are empty?
#' row_na <- row_all_na(df)
#' df[which_(row_na), ]
#'
#' # Removing the empty rows
#' df[which_(row_na, invert = TRUE), ]
#'
#' @rdname is_na
#' @export
is_na <- function(x){
  UseMethod("is_na")
}
#' @rdname is_na
#' @export
is_na.default <- function(x){
  .Call(`_cheapr_cpp_is_na`, x)
}
#' @rdname is_na
#' @export
is_na.POSIXlt <- function(x){
  row_any_na(list_as_df(do.call(recycle, unclass(x)[
    c("sec", "min", "hour", "mday",
      "mon", "year", "wday", "yday")])))
}
#' @rdname is_na
#' @export
is_na.vctrs_rcrd <- function(x){
  row_all_na(list_as_df(x))
}
#' @rdname is_na
#' @export
is_na.data.frame <- function(x){
  row_all_na(x)
}
#' @rdname is_na
#' @export
num_na <- function(x, recursive = TRUE){
  .Call(`_cheapr_cpp_num_na`, x, recursive)
}
#' @rdname is_na
#' @export
which_na <- function(x){
  .Call(`_cheapr_cpp_which_na`, x)
}
#' @rdname is_na
#' @export
which_not_na <- function(x){
  .Call(`_cheapr_cpp_which_not_na`, x)
}
#' @rdname is_na
#' @export
any_na <- function(x, recursive = TRUE){
  .Call(`_cheapr_cpp_any_na`, x, recursive)
}
#' @rdname is_na
#' @export
all_na <- function(x, recursive = TRUE){
  .Call(`_cheapr_cpp_all_na`, x, TRUE, recursive)
}
#' @rdname is_na
#' @export
row_na_counts <- function(x){
  if (!inherits(x, c("matrix", "data.frame"))){
    stop("x must be a matrix or data frame")
  }
  if (is.matrix(x)){
    cpp_matrix_row_na_counts(x)
  } else {
    cpp_row_na_counts(x)
  }
}
#' @rdname is_na
#' @export
col_na_counts <- function(x){
  if (!inherits(x, c("matrix", "data.frame"))){
    stop("x must be a matrix or data frame")
  }
  if (is.matrix(x)){
    cpp_matrix_col_na_counts(x)
  } else {
    cpp_col_na_counts(x)
  }
}

#' @rdname is_na
#' @export
row_all_na <- function(x){
  if (!inherits(x, c("matrix", "data.frame"))){
    stop("x must be a matrix or data frame")
  }
  if (is.matrix(x)){
    cpp_matrix_missing_row(x, 1, TRUE)
  } else {
    cpp_missing_row(x, 1, TRUE)
  }
}
#' @rdname is_na
#' @export
col_all_na <- function(x){
  if (!inherits(x, c("matrix", "data.frame"))){
    stop("x must be a matrix or data frame")
  }
  if (is.matrix(x)){
    cpp_matrix_missing_col(x, 1, TRUE)
  } else {
    cpp_missing_col(x, 1, TRUE)
  }
}
#' @rdname is_na
#' @export
row_any_na <- function(x){
  if (!inherits(x, c("matrix", "data.frame"))){
    stop("x must be a matrix or data frame")
  }
  if (is.matrix(x)){
    cpp_matrix_missing_row(x, 1, FALSE)
  } else {
    cpp_missing_row(x, 1, FALSE)
  }
}
#' @rdname is_na
#' @export
col_any_na <- function(x){
  if (!inherits(x, c("matrix", "data.frame"))){
    stop("x must be a matrix or data frame")
  }
  if (is.matrix(x)){
    cpp_matrix_missing_col(x, 1, FALSE)
  } else {
    cpp_missing_col(x, 1, FALSE)
  }
}

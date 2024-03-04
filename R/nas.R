#' Efficient functions for dealing with missing values.
#'
#' @description
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
#' When `x` is a list, `num_na`, `any_na` and `all_na` will recurse through a
#' potentially nested list for `NA` values. \cr
#'
#' ### Common use-cases
#' To replicate `complete.cases(x)`, use `!row_any_na(x)`. \cr
#' To find rows with any empty values,
#' use `which_(row_any_na(df))`. \cr
#' To find empty rows use `which_(row_all_na(df))`.
#'
#' @param x A vector, matrix or data frame.
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
#' @rdname num_na
#' @export
num_na <- cpp_num_na
#' @rdname num_na
#' @export
which_na <- cpp_which_na
#' @rdname num_na
#' @export
which_not_na <- cpp_which_not_na
#' @rdname num_na
#' @export
any_na <- cpp_any_na
#' @rdname num_na
#' @export
all_na <- function(x){
  cpp_all_na(x, TRUE)
}
#' @rdname num_na
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
#' @rdname num_na
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

#' @rdname num_na
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
#' @rdname num_na
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
#' @rdname num_na
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
#' @rdname num_na
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

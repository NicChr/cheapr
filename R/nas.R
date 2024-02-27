#' Efficient functions for dealing with missing values.
#'
#' @description
#' `num_na(x)` is a faster and more `sum(is.na(x))`. \cr
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
#' All these functions can be parallelised through `options(cheapr.cores)`. \cr
#' When `x` is a data frame, `num_na`, `which_na` and `which_not_na` define a
#' missing value as an empty row with only `NA` values across all columns.
#' To get the number of `NA` values across an entire data frame,
#' use `sum(col_na_counts(data))`. \cr
#' To replicate `complete.cases(x)`, use `!row_any_na(x)`.
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
#' # Number of empty rows
#' num_na(df)
#' # Which rows are empty?
#' which_na(df)
#' df[which_na(df), ]
#'
#' # Removing the empty rows
#' df[which_not_na(df), ]
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
    return(cpp_matrix_col_na_counts(x))
  }
  n_row <- length(attr(x, "row.names"))
  n_col <- length(names(x))
  if (n_row > .Machine$integer.max){
    out <- numeric(n_col)
  } else {
    out <- integer(n_col)
  }
   for (i in seq_len(n_col)){
    out[i] <- num_na(.subset2(x, i))
   }
  # names(out) <- names(x)
  out
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
    n_row <- length(attr(x, "row.names"))
    n_col <- length(names(x))
    out <- logical(n_col)
    for (i in seq_len(n_col)){
      out[i] <- num_na(.subset2(x, i)) == n_row
    }
    # names(out) <- names(x)
    out
  }
  # if (!inherits(x, "data.frame")){
  #   stop("x must be a data frame")
  # }
  # n_row <- length(attr(x, "row.names"))
  # n_col <- length(names(x))
  # out <- logical(n_col)
  # for (i in seq_len(n_col)){
  #   out[i] <- num_na(.subset2(x, i)) == n_row
  # }
  # # names(out) <- names(x)
  # out
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
    n_row <- length(attr(x, "row.names"))
    n_col <- length(names(x))
    out <- logical(n_col)
    for (i in seq_len(n_col)){
      out[i] <- num_na(.subset2(x, i)) > 0
    }
    out
  }
}

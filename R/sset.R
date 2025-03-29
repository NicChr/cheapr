#' Cheaper subset
#'
#' @description
#' Cheaper alternative to `[` that consistently subsets data frame
#' rows, always returning a data frame. There are explicit methods for
#' enhanced data frames like tibbles, data.tables and sf.
#'
#' @param x Vector or data frame.
#' @param i A logical or vector of indices. \cr
#' @param j Column indices, names or logical vector.
#' @param ... Further parameters passed to `[`.
#'
#' @returns
#' A new vector, data frame, list, matrix or other R object.
#'
#' @details
#' `sset` is an S3 generic.
#' You can either write methods for `sset` or `[`. \cr
#' `sset` will fall back on using `[` when no suitable method is found.
#'
#' To get into more detail, using `sset()` on a data frame, a new
#' list is always allocated through `new_list()`.
#'
#' ### Difference to base R
#'
#' When `i` is a logical vector, it is passed directly to `which_()`. \cr
#' This means that `NA` values are ignored and this also means that `i`
#' is not recycled, so it is good practice to make sure the logical vector
#' matches the length of x. To return `NA` values, use `sset(x, NA_integer_)`.
#'
#' ### ALTREP range subsetting
#'
#' When `i` is an ALTREP compact sequence which can be commonly created
#' using e.g. `1:10` or using `seq_len`, `seq_along` and `seq.int`,
#' `sset` internally uses a range-based subsetting method which is faster and doesn't
#' allocate `i` into memory.
#'
#' @examples
#' library(cheapr)
#' library(bench)
#'
#' # Selecting columns
#' sset(airquality, j = "Temp")
#' sset(airquality, j = 1:2)
#'
#' # Selecting rows
#' sset(iris, 1:5)
#'
#' # Rows and columns
#' sset(iris, 1:5, 1:5)
#' sset(iris, iris$Sepal.Length > 7, c("Species", "Sepal.Length"))
#'
#' # Comparison against base
#' x <- rnorm(10^4)
#'
#' mark(x[1:10^3], sset(x, 1:10^3))
#' mark(x[x > 0], sset(x, x > 0))
#'
#' df <- data.frame(x = x)
#'
#' mark(df[df$x > 0, , drop = FALSE],
#'      sset(df, df$x > 0),
#'      check = FALSE) # Row names are different
#'
#' @rdname sset
#' @export
sset <- function(x, ...){
  UseMethod("sset")
}
#' @export
sset.default <- function(x, i, ...){
  if (cpp_is_simple_vec(x) && n_dots(...) == 0){
    .Call(`_cheapr_cpp_sset`, x, if (missing(i)) seq_len(vector_length(x)) else i, TRUE)
  } else {
    if (!missing(i) && is.logical(i)){
      check_length(i, length(x))
      i <- which_(i)
    }
    x[i, ...]
  }
}
#' @rdname sset
#' @export
sset.data.frame <- function(x, i = NULL, j = NULL, ...){
  .Call(`_cheapr_cpp_df_subset`, x, i, j)
}
#' @rdname sset
#' @export
sset.POSIXlt <- function(x, i = NULL, j = NULL, ...){
  missingj <- is.null(j)
  out <- fill_posixlt(x, classed = FALSE)
  if (missingj){
    j <- seq_along(out)
  }
  out <- sset_df(list_as_df(out), i , j)
  if (missingj){
    set_attr(out, "class", class(x))
    set_rm_attr(out, "row.names")
  }
  set_attr(out, "tzone", attr(x, "tzone"))
  if (posixlt_is_balanced(x)){
    set_attr(out, "balanced", TRUE)
  } else {
    set_attr(out, "balanced", NA)
  }
  out
}
#' @rdname sset
#' @export
sset.sf <- function(x, i = NULL, j = NULL, ...){
  out <- sset_df(x, i, j)
  cpp_reconstruct(
    out, x, c("names", "row.names"), cpp_setdiff(names(attributes(x)), c("names", "row.names"))
  )
}
#' @export
sset.vctrs_rcrd <- function(x, i = NULL, ...){
  out <- sset_row(list_as_df(x), i)
  attributes(out) <- attributes(x)
  out
}

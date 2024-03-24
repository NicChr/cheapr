#' Cheaper subset
#'
#' @description
#' Cheaper alternative to `[` that consistently subsets data frame
#' rows, always returning a data frame. There are explicit methods for
#' enhanced data frames like tibbles, data.tables and sf.
#'
#' @param x Vector or data frame.
#' @param i A logical or vector of indices. \cr
#' The default is 0 which differs to `[`.
#' @param j Column indices, names or logical vector.
#' @param ... Further parameters passed to `[`.
#'
#' @details
#' `sset` is an S3 generic.
#' You can either write methods for `sset` or `[`. \cr
#' `sset` will fall back on using `[` when no suitable method is found.
#'
#' @examples
#' library(cheapr)
#' library(bench)
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
sset <- function(x, i = 0, ...){
  UseMethod("sset")
}
#' @export
sset.default <- function(x, i = 0, ...){
  if (is.logical(i)){
    check_length(i, length(x))
    i <- which_(i)
  }
  x[i, ...]
}
#' @rdname sset
#' @export
sset.data.frame <- function(x, i = 0, j = seq_along(x), ...){
  df_subset(x, i, j)
}
#' @rdname sset
#' @export
sset.tbl_df <- function(x, i = 0, j = seq_along(x), ...){
  out <- df_subset(x, i, j)
  class(out) <- c("tbl_df", "tbl", "data.frame")
  out
}
#' @rdname sset
#' @export
sset.POSIXlt <- function(x, i = 0, ...){
  out <- df_subset(list_as_df(x), i)
  cpp_set_copy_attributes(out, x, names(attributes(x)))
}
#' @rdname sset
#' @export
sset.data.table <- function(x, i = 0, j = seq_along(x), ...){
  # collapse::qDT(df_subset(x, i, j))
  out <- df_subset(x, i, j)
  cpp_set_copy_attributes(
    out, x, c("class", ".internal.selfref")
  )
}
#' @rdname sset
#' @export
sset.sf <- function(x, i = 0, j = seq_along(x), ...){
  out <- df_subset(x, i, j)
  source_nms <- names(attributes(x))
  invisible(
    cpp_set_copy_attributes(out, x, setdiff_(source_nms, c("names", "row.names", "class")))
  )
  class(out) <- class(x)
  out
}

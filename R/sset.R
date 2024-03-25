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
#' @details
#' `sset` is an S3 generic.
#' You can either write methods for `sset` or `[`. \cr
#' `sset` will fall back on using `[` when no suitable method is found.
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
  if (!missing(i) && is.logical(i)){
    check_length(i, length(x))
    i <- which_(i)
  }
  x[i, ...]
}
#' @rdname sset
#' @export
sset.data.frame <- function(x, i, j = seq_along(x), ...){
  df_subset(x, i, j)
}
#' @rdname sset
#' @export
sset.tbl_df <- function(x, i, j = seq_along(x), ...){
  out <- df_subset(x, i, j)
  class(out) <- c("tbl_df", "tbl", "data.frame")
  out
}
#' @rdname sset
#' @export
sset.POSIXlt <- function(x, i, ...){
  out <- df_subset(list_as_df(x), i)
  cpp_set_copy_attributes(
    cpp_set_rm_attributes(out), x, names(attributes(x))
  )
}
#' @rdname sset
#' @export
sset.data.table <- function(x, i, j = seq_along(x), ...){
  # This is to ensure that a copy is made basically
  # More efficient to use data.table::copy()
  if (missing(i)){
    i <- seq_len(nrow(x))
  }
  out <- df_subset(x, i, j)
  cpp_set_copy_attributes(
    out, x, c("class", ".internal.selfref")
  )
}
#' @rdname sset
#' @export
sset.sf <- function(x, i, j = seq_along(x), ...){
  out <- df_subset(x, i, j)
  source_nms <- names(attributes(x))
  invisible(
    cpp_set_copy_attributes(out, x, setdiff_(source_nms, c("names", "row.names", "class")))
  )
  class(out) <- class(x)
  out
}
df_select <- function(x, j){
  if (is.logical(j)){
    check_length(j, length(x))
    j <- which_(j)
  }
  if (is.character(j)){
    j <- collapse::fmatch(j, names(x), overid = 2L)
  }
  attrs <- attributes(x)
  out <- cpp_list_rm_null(unclass(x)[j])
  attrs[["names"]] <- attr(out, "names")
  attrs[["row.names"]] <- .row_names_info(x, type = 0L)
  attributes(out) <- attrs
  out
}

# Efficient data frame subset
# With the exception of which_() this is surprisingly all base R...
# It relies on sset which falls back on `[` when no method is found.
df_subset <- function(x, i, j = seq_along(x)){
  missingi <- missing(i)
  nrows <- length(attr(x, "row.names"))
  if (!missingi && is.logical(i)){
    check_length(i, nrows)
    i <- which_(i)
  }

  ### Subset columns

  out <- df_select(x, j)

  ### Subset rows

  if (!missingi){
    if (length(out) == 0){
      attr(out, "row.names") <- .set_row_names(
        length(
          seq_len(nrows)[i]
        )
      )
    } else {
      out <- list_as_df(
        lapply(out, sset, i)
      )
    }
  }
  out
}

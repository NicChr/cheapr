#' Cheaper subset `sset()`
#'
#' @description
#'
#' `sset()` is a cheaper alternative to `[`.
#'
#' It consistently subsets data frame rows for any data frame class including
#' tibble and data.table.
#'
#' @param x Vector or data frame.
#' @param i A logical vector or integer vector of locations.
#' @param j Column indices, names or logical vector.
#' @param ... Further parameters passed to `[`.
#'
#' @returns
#' A new vector, data frame, list, matrix or other R object.
#'
#' @details
#'
#' ### S3 dispatching
#' `sset` will internally dispatch the correct method and
#' will call `[` if it can't find an appropriate method. This means
#' one can define their own `[` method for custom S3 objects.
#'
#' To speed up subsetting for common objects likes Dates and `POSIXlt`
#' an internal generic function is used which overwrites the `[` method for
#' that common object. This is why subsetting `POSIXlt` is much faster with
#' `sset` an internal method has been defined. For more details see the code
#' for `cheapr:::cheapr_sset`.
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
#' `sset` internally uses a range-based subsetting method
#' which is faster and doesn't allocate `i` into memory.
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
sset <- function(x, i = NULL, j = NULL, ...){
  .Call(`_cheapr_cpp_sset2`, x, i, j, TRUE, list(...))
}
cheapr_sset <- function(x, ...){
  UseMethod("cheapr_sset")
}
#' @export
cheapr_sset.default <- function(x, i = NULL, j = NULL, ...){
  if (is.logical(i)){
    check_length(i, length(x))
    i <- which_(i)
  }
  if (is.null(i) && is.null(j)){
    x[...]
  } else if (is.null(j)){
    x[i = i, ...]
  } else if (is.null(i)){
    x[j = j, ...]
    } else {
    x[i = i, j = j, ...]
  }
}
#' @export
cheapr_sset.data.frame <- function(x, i = NULL, j = NULL, ...){
  sset_df(x, i, j, ...)
}
#' @export
cheapr_sset.POSIXlt <- function(x, i = NULL, j = NULL, ...){
  missingj <- is.null(j)
  out <- fill_posixlt(x, classed = FALSE)
  if (missingj){
    j <- seq_along(out)
  }
  out <- sset(list_as_df(out), i, j)
  if (missingj){
    attrs_add(out, class = class(x), row.names = NULL, .set = TRUE)
  }
  attrs_add(out, tzone = attr(x, "tzone"), .set = TRUE)
  if (posixlt_is_balanced(x)){
    attrs_add(out, balanced = TRUE, .set = TRUE)
  } else {
    attrs_add(out, balanced = NA, .set = TRUE)
  }
  out
}
# Have to define sf method because otherwise `cheapr_sset.data.frame` is called
# which can't handle 'sf_column' correctly
#' @export
cheapr_sset.sf <- function(x, i = NULL, j = NULL, ...){
  cheapr_sset.default(x, i, j, ...)
}
#' @export
cheapr_sset.vctrs_rcrd <- function(x, i = NULL, ...){
  x |>
    list_as_df() |>
    sset(i) |>
    attrs_clear(.set = TRUE) |>
    attrs_modify(.args = shallow_copy(attributes(x)), .set = TRUE)
}

#' @export
cheapr_sset.integer64 <- function(x, i = NULL, ...){
  cpp_sset_int64(x, i) |>
    attrs_modify(class = "integer64", .set = TRUE)
}

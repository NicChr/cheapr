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
#' To get into more detail, using `sset()` on a data frame, a new
#' list is always allocated through `cheapr:::cpp_new_list()`.
#'
#' ### Difference to base R
#'
#' When `i` is a logical vector, it is passed directly to `which_()`. \cr
#' This means that `NA` values are ignored and this also means that `i`
#' is not recycled, so it is good practice to make sure the logical vector
#' matches the length of x. To return `NA` values, use `x[NA_integer_]`.
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
sset.POSIXlt <- function(x, i, j, ...){
  missingi <- missing(i)
  missingj <- missing(j)
  if (n_unique(lengths_(unclass(x))) > 1){
    out <- balancePOSIXlt(x, fill.only = FALSE, classed = FALSE)
  } else {
    out <- unclass(x)
  }
  if (missingj){
    j <- seq_along(out)
  }
  out <- df_subset(list_as_df(out), i, j)
  cpp_set_rm_attr(out, "row.names")
  if (missingj){
    cpp_set_add_attr(out, "class", class(x))
  }
  cpp_set_add_attr(out, "tzone", attr(x, "tzone"))
  cpp_set_add_attr(out, "balanced", TRUE)
}
#' @rdname sset
#' @export
sset.data.table <- function(x, i, j = seq_along(x), ...){
  out <- df_subset(x, i, j)
  cpp_set_add_attributes(out, list(class = class(x),
                               .internal.selfref = attributes(x)[[".internal.selfref"]]),
                     add = TRUE)
  dt_alloc <- tryCatch(get("setalloccol",
                           asNamespace("data.table"),
                           inherits = FALSE),
                       error = function(e) return(".r.error"))
  # Reserve sufficient space as data.table::truelength(out) at this point is 0
  if (is.character(dt_alloc) && length(dt_alloc) == 1 && dt_alloc == ".r.error"){
    out <- collapse::qDT(out)
  } else {
    dt_alloc(out, n = getOption("datatable.alloccol", 1024L))
  }
  out
}
#' @rdname sset
#' @export
sset.sf <- function(x, i, j = seq_along(x), ...){
  out <- df_subset(x, i, j)
  source_attrs <- attributes(x)
  source_nms <- names(source_attrs)
  attrs_to_keep <- source_attrs[setdiff_(source_nms, c("names", "row.names"))]
  cpp_set_add_attributes(out, attrs_to_keep, add = TRUE)
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
  cpp_set_add_attributes(out, attrs, add = FALSE)
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
    i <- as.integer(i)
    if (length(out) == 0){
      attr(out, "row.names") <- .set_row_names(
        length(
          seq_len(nrows)[i]
        )
      )
    } else {
      out <- cpp_df_sset(out, i)
      # out <- list_as_df(
      #   lapply(out, sset, i)
      # )
    }
  }
  out
}
# Turn negative indices to positives
neg_indices_to_pos <- function(exclude, n){
  if (n == 0){
    integer()
  } else {
    which_not_in(seq.int(from = -1L, to = -n, by = -1L), exclude)
  }
  # out <- which_not_in(seq_len(n), cpp_set_reverse_sign(exclude))
  # exclude <- cpp_set_reverse_sign(exclude)
  # out
  # which_not_in(seq_len(n), abs(exclude))
}
clean_indices <- function(indices, n){
  valid <- cpp_index_is_valid(indices, n)
  if (count_val(valid, TRUE) == length(indices)){
    indices
  } else {
    indices[which_(valid)]
  }
}
# Efficient function to remove a value from a vector
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

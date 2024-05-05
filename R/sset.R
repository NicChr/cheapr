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
  if (!missing(i) && is.logical(i)){
    check_length(i, length(x))
    i <- which_(i)
  }
  # The below line will handle a special but common
  # case of subsetting with a fairly large altrep int sequence
  # For non-classed objects
  if (!is.object(x) && !missing(i) &&
      is_compact_seq(i) && n_dots(...) == 0){
    int_seq_data <- compact_seq_data(i)
    from <- int_seq_data[[1]]
    to <- int_seq_data[[2]]
    by <- int_seq_data[[3]]
    out <- cpp_sset_range(x, from, to, by)
    if (!is.null(names(x))){
      names(out) <- cpp_sset_range(names(x), from, to, by)
    }
    out
  } else {
    x[i, ...]
  }
}
#' @rdname sset
#' @export
sset.Date <- function(x, i, ...){
  if (!missing(i) && is.logical(i)){
    check_length(i, length(x))
    i <- which_(i)
  }
  if (!missing(i) &&
      is_compact_seq(i) && n_dots(...) == 0){
    int_seq_data <- compact_seq_data(i)
    from <- int_seq_data[[1]]
    to <- int_seq_data[[2]]
    by <- int_seq_data[[3]]
    out <- cpp_sset_range(x, from, to, by)
    if (!is.null(names(x))){
      set_attr(out, "names", cpp_sset_range(names(x), from, to, by))
    }
    set_attr(out, "class", oldClass(x))

  } else {
    x[i, ...]
  }
}
#' @rdname sset
#' @export
sset.POSIXct <- function(x, i, ...){
  if (!missing(i) && is.logical(i)){
    check_length(i, length(x))
    i <- which_(i)
  }
  if (!missing(i) &&
      is_compact_seq(i) && n_dots(...) == 0){
    int_seq_data <- compact_seq_data(i)
    from <- int_seq_data[[1]]
    to <- int_seq_data[[2]]
    by <- int_seq_data[[3]]
    out <- cpp_sset_range(x, from, to, by)
    if (!is.null(names(x))){
      set_attr(out, "names", cpp_sset_range(names(x), from, to, by))
    }
    set_attr(out, "tzone", attr(x, "tzone"))
    set_attr(out, "class", oldClass(x))
  } else {
    x[i, ...]
  }
}
#' @rdname sset
#' @export
sset.factor <- function(x, i, ...){
  if (!missing(i) && is.logical(i)){
    check_length(i, length(x))
    i <- which_(i)
  }
  if (!missing(i) &&
      is_compact_seq(i) && n_dots(...) == 0){
    int_seq_data <- compact_seq_data(i)
    from <- int_seq_data[[1]]
    to <- int_seq_data[[2]]
    by <- int_seq_data[[3]]
    out <- cpp_sset_range(x, from, to, by)
    if (!is.null(names(x))){
      set_attr(out, "names", cpp_sset_range(names(x), from, to, by))
    }
    set_attr(out, "levels", attr(x, "levels"))
    set_attr(out, "class", oldClass(x))
  } else {
    x[i, ...]
  }
}
#' @rdname sset
#' @export
sset.data.frame <- function(x, i, j, ...){
  df_subset(x, i, j)
}
#' @rdname sset
#' @export
sset.tbl_df <- function(x, i, j, ...){
  out <- df_subset(x, i, j)
  class(out) <- c("tbl_df", "tbl", "data.frame")
  out
}
#' @rdname sset
#' @export
sset.POSIXlt <- function(x, i, j, ...){
  missingi <- missing(i)
  missingj <- missing(j)
  out <- fill_posixlt(x, classed = FALSE)
  if (missingj){
    j <- seq_along(out)
  }
  out <- df_subset(list_as_df(out), i, j)
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
sset.data.table <- function(x, i, j, ...){
  out <- df_subset(x, i, j)
  set_attrs(out, list(
    class = class(x),
    .internal.selfref = attributes(x)[[".internal.selfref"]]
  ), add = TRUE)
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
sset.sf <- function(x, i, j, ...){
  out <- df_subset(x, i, j)
  source_attrs <- attributes(x)
  source_nms <- names(source_attrs)
  attrs_to_keep <- source_attrs[setdiff_(source_nms, c("names", "row.names"))]
  set_attrs(out, attrs_to_keep, add = TRUE)
}
df_select <- function(x, j){
  j_exists <- !missing(j)
  if (j_exists && is.logical(j)){
    check_length(j, length(x))
    j <- which_(j)
  }
  if (j_exists && is.character(j)){
    j <- collapse::fmatch(j, names(x), overid = 2L)
  }
  attrs <- attributes(x)
  if (j_exists){
    out <- cpp_list_rm_null(.subset(x, j), always_shallow_copy = FALSE)
  } else {
    out <- cpp_list_rm_null(x)
  }
  # Neater but not as efficient for dfs with many cols
  # out <- cpp_list_rm_null(unclass(x)[j])
  attrs[["names"]] <- attr(out, "names")
  attrs[["row.names"]] <- .row_names_info(x, type = 0L)
  set_attrs(out, attrs, add = FALSE)
}

# Efficient data frame subset
# It relies on sset which falls back on `[` when no method is found.
df_subset <- function(x, i, j){
  missingi <- missing(i)
  missingj <- missing(j)
  nrows <- length(attr(x, "row.names"))
  if (!missingi && is.logical(i)){
    check_length(i, nrows)
    i <- which_(i)
  }

  ### Subset columns
  # If j arg is missing, we want to still create a shallow copy
  # Which we do through df_select()
  if (!missingj || (missingi && missingj)){
    out <- df_select(x, j)
  } else {
    out <- x
  }
  ### Subset rows
  if (!missingi){
    out <- cpp_sset_df(out, as.integer(i))
  }
  out
}
# Turn negative indices to positives
neg_indices_to_pos <- function(exclude, n){
  if (n == 0){
    integer()
  } else {
    which_not_in(
      seq.int(from = -1L, to = -as.integer(n), by = -1L),
      as.integer(exclude)
    )
  }
}

#' Cheap data frame utilities
#'
#' @param ... Key-value pairs.
#' @param .nrows `[integer(1)]` - (Optional) number of rows. \cr
#' Commonly used to initialise a 0-column data frame with rows.
#' @param .recycle `[logical(1)]` - Should arguments be recycled?
#' Default is `TRUE`.
#' @param .name_repair `[logical(1)]` - Should duplicate and empty names
#' repaired and made unique? Default is `TRUE`.
#' @param x An object to coerce to a `data.frame` or a character vector
#' for `unique_name_repair()`.
#' @param dup_sep `[character(1)]` A separator to use between
#' duplicate column names and their locations. Default is `'_'`
#' @param empty_sep `[character(1)]` A separator to use between the empty
#' column names and their locations. Default is `'col_'`
#' @param .args An alternative to `...` for easier programming with lists.
#' @param cols A list of values to add or modify data frame `x`.
#'
#' @returns
#' A `data.frame`. \cr
#' `name_repair` takes a character vector and returns unique strings by
#' appending duplicate string locations to the duplicates.
#' This is mostly used to create unique col names.
#'
#' @details
#' `fast_df()` is a very fast bare-bones version of `new_df()` that
#' performs no checks and no recycling or name tidying, making it appropriate
#' for very tight loops.
#'
#' @rdname data_frame
#' @export
new_df <- function(..., .nrows = NULL, .recycle = TRUE, .name_repair = TRUE, .args = NULL){
  .Call(`_cheapr_cpp_new_df`, .Call(`_cheapr_cpp_list_args`, list(...), .args),
        .nrows, .recycle, .name_repair)
}
#' @rdname data_frame
#' @export
as_df <- function(x){
  if (inherits(x, "data.frame")){
    return(cpp_new_df(x, length(attr(x, "row.names")), FALSE, FALSE))
  } else if (is.null(x) || (is.atomic(x) && length(dim(x)) < 2)){
    out <- list_drop_null(list(name = names(x), value = x))
    attr(out, "row.names") <- .set_row_names(NROW(x))
    class(out) <- "data.frame"
    return(out)
  } else {
    # Plain list
    if (!is.object(x) && is.list(x)){
      return(new_df(.args = x, .recycle = TRUE, .name_repair = TRUE))
    } else {
      out <- as.data.frame(x, stringsAsFactors = FALSE)
      names(out) <- str_coalesce(names(out), paste0("col_", seq_along(out)))
      return(out)
    }
  }
}

#' @rdname data_frame
#' @export
fast_df <- function(..., .args = NULL){
  .Call(`_cheapr_cpp_list_as_df`, .Call(`_cheapr_cpp_list_args`, list(...), .args))
}

#' @rdname data_frame
#' @export
df_modify <- cpp_df_assign_cols

#' @rdname data_frame
#' @export
list_as_df <- cpp_list_as_df

#' @rdname data_frame
#' @export
name_repair <- function(x, dup_sep = "_", empty_sep = "col_"){
  cpp_name_repair(x, dup_sep, empty_sep)
}
#' @rdname data_frame
#' @export
unique_name_repair <- name_repair
#' @rdname data_frame
#' @export
col_c <- function(..., .recycle = TRUE, .name_repair = TRUE, .args = NULL){
  .Call(
    `_cheapr_cpp_df_col_c`, .Call(`_cheapr_cpp_list_args`, list(...), .args),
    .recycle, .name_repair
  )
}
#' @rdname data_frame
#' @export
row_c <- cheapr_c

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
#' @param .args An alternative to `...` so you can supply arguments directly
#' in a list. \cr
#' This is equivalent to `do.call(f, .args)` but much more efficient.
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
as_df <- cpp_as_df

#' @rdname data_frame
#' @export
fast_df <- function(..., .args = NULL){
  .Call(`_cheapr_cpp_new_df`, .Call(`_cheapr_cpp_list_args`, list(...), .args), NULL, FALSE, FALSE)
}

#' @rdname data_frame
#' @export
df_modify <- cpp_df_assign_cols

#' @rdname data_frame
#' @export
list_as_df <- function(x){
  .Call(`_cheapr_cpp_new_df`, x, NULL, FALSE, FALSE)
}

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

#' Fast functions for data frame subsetting
#'
#' @description
#' These functions are for developers that need minimal overhead when
#' filtering on rows and/or cols.
#'
#' @param x A `data.frame`.
#' @param i Rows - If `NULL` all rows are returned.
#' @param j Cols - If `NULL` all cols are returned.
#' @param ... Unused.
#'
#' @returns
#' A data frame subsetted on rows `i` and cols `j`.
#'
#' @details
#' If you are unsure which functions to use then it is recommended to use
#' `sset()`. These low-overhead helpers do not work well with data.tables
#' but should work well with basic data frames and basic tibbles.
#' The only real difference between `sset_df` and `sset_row`/`sset_col` is that
#' `sset_df` attempts to return a similar type of data frame as the input,
#' whereas `sset_row` and `sset_col` always return a plain data frame.
#'
#' @rdname sset_df
#' @export
sset_df <- function(x, i = NULL, j = NULL, ...){
  .Call(`_cheapr_cpp_df_subset`, x, i, j, TRUE)
}
#' @rdname sset_df
#' @export
sset_row <- function(x, i = NULL){
  .Call(`_cheapr_cpp_df_slice`, x, i, TRUE)
}
#' @rdname sset_df
#' @export
sset_col <- function(x, j = NULL){
  .Call(`_cheapr_cpp_df_select`, x, j)
}

# Keep this for fastplyr otherwise it breaks dependency
df_select <- function(x, j = NULL){
  cpp_reconstruct(
    sset_col(x, missing(j) %!||% j), x,
    "names", val_rm(names(attributes(x)), "names"),
    TRUE
  )
}

# Kept for reverse compatibility reasons
cpp_sset_df <- sset_row

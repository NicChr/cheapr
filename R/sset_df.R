#' Fast functions for data frame subsetting
#'
#' @description
#' These functions are for developers that need minimal overhead when
#' filtering on rows and/or cols.
#'
#' @param x A `data.frame`.
#' @param i Rows - If `NULL` all rows are returned.
#' @param j Cols - If `NULL` all cols are returned.
#' @param keep_attrs Should all attributes (except for `names` and `row.names`)
#' be kept as is? The default is `FALSE` which returns a plain data frame.
#'
#' @returns
#' A data frame subsetted on rows `i` and cols `j`.
#'
#' @details
#' If you are unsure which functions to use then it is recommended to use
#' `sset()`. These low-overhead helpers do not work well with data.tables
#' but should work well with basic data frames and basic tibbles.
#'
#' @rdname sset_df
#' @export
sset_df <- function(x, i = NULL, j = NULL, keep_attrs = FALSE){
  .Call(`_cheapr_cpp_df_subset`, x, i, j, keep_attrs)
}
#' @rdname sset_df
#' @export
sset_row <- function(x, i = NULL){
  .Call(`_cheapr_cpp_df_slice`, x, i)
}
#' @rdname sset_df
#' @export
sset_col <- function(x, j = NULL){
  .Call(`_cheapr_cpp_df_select`, x, j)
}

# Keep this for fastplyr otherwise it breaks dependency
df_select <- function(x, j = NULL){
  sset_df(x, i = NULL, j = missing(j) %!||% j, keep_attrs = TRUE)
}

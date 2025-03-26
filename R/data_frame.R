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
#' column names and their locations. Default is `'_'`
#'
#' @returns
#' A `data.frame`. \cr
#' `unique_name_repair` takes a character vector and returns unique strings by
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
new_df <- function(..., .nrows = NULL, .recycle = TRUE, .name_repair = TRUE){
  .Call(`_cheapr_cpp_new_df`, list(...), .nrows, .recycle, .name_repair)

}
#' @rdname data_frame
#' @export
as_df <- function(x){
  if (inherits(x, "data.frame")){
    out <- x
    class(out) <- "data.frame"
  } else if (is.null(x) || (is.atomic(x) && length(dim(x)) < 2)){
    out <- cpp_list_rm_null(list(name = names(x), value = x))
    attr(out, "row.names") <- .set_row_names(NROW(x))
    class(out) <- "data.frame"
  } else {
    # Plain list
    if (!is.object(x) && is.list(x)){
      out <- list_as_df(cpp_recycle(x, NULL))
    } else {
      out <- as.data.frame(x, stringsAsFactors = FALSE)
    }
    if (is.null(names(out))){
      names(out) <- paste0("col_", seq_along(out))
    }
    non_empty <- nzchar(names(out))
    if (!all(non_empty)){
      empty <- which_(non_empty, invert = TRUE)
      names(out)[empty] <- paste0("col_", empty)
    }
  }
  out
}

#' @rdname data_frame
#' @export
fast_df <- function(...){
  .Call(`_cheapr_cpp_list_as_df`, list(...))
}

#' @rdname data_frame
#' @export
unique_name_repair <- function(x, dup_sep = "_", empty_sep = "_"){
  cpp_name_repair(x, dup_sep, empty_sep)
}

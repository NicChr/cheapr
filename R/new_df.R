#' Fast data frame constructors
#'
#' @param ... Key-value pairs.
#' @param .nrows `[integer(1)]` - (Optional) number of rows. \cr
#' Commonly used to initialise a 0-column data frame with rows.
#' @param .recycle `[logical(1)]` - Should arguments be recycled?
#' Default is `FALSE`.
#' @param .name_repair `[logical(1)]` - Should duplicate names be made unique?
#' Default is `FALSE`.
#' @param x An object to coerce to a `data.frame`.
#'
#' @returns A `data.frame`
#'
#' @details
#' `fast_df()` is a very fast bare-bones version of `new_df()` that
#' performs no checks and no recycling or name tidying.
#' All variables must be named and of equal length.
#'
#' @rdname new_df
#' @export
new_df <- function(..., .nrows = NULL, .recycle = FALSE, .name_repair = FALSE){

  out <- named_list(..., .keep_null = FALSE)

  # Recycle
  if (.recycle){
    out <- cpp_recycle(out, length = .nrows)
  }

  if (is.null(.nrows)){
    if (length(out) == 0L){
      row_names <- integer()
    } else {
      N <- NROW(.subset2(out, 1L))
      row_names <- c(NA_integer_, -N)
    }
  } else {
    row_names <- .set_row_names(.nrows)
  }

  out_names <- as.character(attr(out, "names", TRUE))

  if (.name_repair){
    out_names <- unique_name_repair(out_names)
  }

  attr(out, "names") <- out_names
  attr(out, "row.names") <- row_names
  class(out) <- "data.frame"
  out
}
#' @rdname new_df
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
      out <- list_as_df(do.call(recycle, as.list(x)))
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
unique_name_repair <- function(x, .sep = "..."){
  cpp_name_repair(x, .sep)
}

#' @rdname new_df
#' @export
fast_df <- function(...){
  .Call(`_cheapr_cpp_list_as_df`, list(...))
}

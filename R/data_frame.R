#' Cheap data frame utilities
#'
#' @param ... Key-value pairs.
#' @param .nrows `[integer(1)]` - (Optional) number of rows. \cr
#' Commonly used to initialise a 0-column data frame with rows.
#' @param .recycle `[logical(1)]` - Should arguments be recycled?
#' Default is `FALSE`.
#' @param .name_repair `[logical(1)]` - Should duplicate names be made unique?
#' Default is `FALSE`.
#' @param x An object to coerce to a `data.frame` or a character vector
#' for `unique_name_repair()`.
#' @param dup_sep `[character(1)]` A separator to use between
#' duplicate column names and their locations. Default is `'_'`
#' @param empty_sep `[character(1)]` A separator to use between the empty
#' column names and their locations. Default is `'_'`
#'
#' @returns A `data.frame`. \cr
#' `unique_name_repair` takes a character vector and returns unique strings by
#' appending duplicate string locations to the duplicates.
#' This is mostly used to create unique col names.
#'
#' @details
#' `fast_df()` is a very fast bare-bones version of `new_df()` that
#' performs no checks and no recycling or name tidying.
#' All variables must be named and of equal length.
#'
#' @rdname data_frame
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

# df_reconstruct <- function(data, template, data_attrs_keep = "", template_attrs_keep = "class"){
#   cpp_df_reconstruct(data, template, data_attrs_keep, template_attrs_keep)
# }

# df_reconstruct <- function(data, template, ...){
#   UseMethod("df_reconstruct", template)
# }
# df_reconstruct.default <- function(data, template, keep_attrs = FALSE, ...){
#   cpp_df_reconstruct(data, template, keep_attrs)
# }

# df_reconstruct.data.table <- function(data, template, keep_attrs = FALSE, copy_all_attrs = FALSE, ...){
#
#   at <- attributes(template)
#   row_names <- .row_names_info(data, type = 0L)
#   out <- collapse::qDT(shallow_copy(data), keep.attr = keep_attrs)
#
#   selfref <- attr(out, ".internal.selfref", TRUE)
#   cpp_set_add_attr(out, "row.names", row_names)
#
#   if (copy_all_attrs){
#     for (a in cpp_setdiff(
#       names(at), c("row.names", "names", ".internal.selfref")
#     )){
#       cpp_set_add_attr(out, a, at[[a]])
#     }
#   } else {
#     cpp_set_add_attr(out, "class", at[["class"]])
#     cpp_set_add_attr(out, ".internal.selfref", selfref)
#   }
#   out
#
#   # at <- attributes(template)
#   # row_names <- .row_names_info(data, type = 0L)
#   # out <- collapse::qDT(shallow_copy(data), keep.attr = keep_attrs)
#   # selfref <- attr(out, ".internal.selfref", TRUE)
#   #
#   # if (!keep_attrs){
#   #   cpp_set_rm_attributes(out)
#   # }
#   # cpp_set_add_attr(out, "names", attr(data, "names", TRUE))
#   # cpp_set_add_attr(out, "row.names", row_names)
#   # if (copy_all_attrs){
#   #   for (a in cpp_setdiff(
#   #     names(at), c("row.names", "names", ".internal.selfref")
#   #   )){
#   #     cpp_set_add_attr(out, a, at[[a]])
#   #   }
#   # } else {
#   #   cpp_set_add_attr(out, "class", at[["class"]])
#   #   cpp_set_add_attr(out, ".internal.selfref", selfref)
#   # }
#   # out
# }

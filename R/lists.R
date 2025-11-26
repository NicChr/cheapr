#' List utilities
#'
#' @description
#' Functions to help work with lists.
#'
#' @param x A list.
#' @param length Length of list.
#' @param default Default value for each list element.
#' @param names Should names of list elements be added? Default is `FALSE`.
#' @param values A named list
#' @param ... Objects to combine into a list.
#' @param .args An alternative to `...` so you can supply arguments directly
#' in a list. \cr
#' This is equivalent to `do.call(f, .args)` but much more efficient.
#'
#' @returns
#' `list_lengths()` returns the list lengths. \cr
#' `unlisted_length()` is a fast alternative to `length(unlist(x))`. \cr
#' `new_list()` is like `vector("list", length)` but also allows you to specify
#' a default value for each list element. This can be useful for
#' initialising with a catch-all value so that when you unlist you're guaranteed
#' a list of length >= to the specified length.
#'
#' `list_assign()` is vectorised version of `[[<-` that
#' concatenates `values` to `x` or modifies `x` where the
#' names match. Can be useful for modifying data frame variables.
#'
#' `list_combine()` combines each element of a set of lists into a single list.
#' If an element is not a list, it is treated as a length-one list.
#' This happens to be very useful for combining data frame cols.
#'
#'
#' `list_drop_null()` removes `NULL` list elements very quickly.
#'
#' @examples
#' library(cheapr)
#' l <- list(1:10,
#'           NULL,
#'           list(integer(), NA_integer_, 2:10))
#'
#' lengths_(l) # Faster lengths()
#' unlisted_length(l) # length of vector if we unlist
#' paste0("length: ", length(print(unlist(l))))
#'
#' unlisted_length(l) - na_count(l) # Number of non-NA elements
#'
#' # We can create and initialise a new list with a default value
#' l <- new_list(20, 0L)
#' l[1:5]
#' # This works well with vctrs_list_of objects
#'
#' @export
#' @rdname lists
list_lengths <- function(x, names = FALSE){
  .Call(`_cheapr_cpp_lengths`, x, names)
}
#' @export
#' @rdname lists
lengths_ <- list_lengths
#' @export
#' @rdname lists
unlisted_length <- cpp_unnested_length
#' @export
#' @rdname lists
new_list <- function(length = 0L, default = NULL){
  .Call(`_cheapr_cpp_new_list`, length, default)
}
#' @export
#' @rdname lists
list_assign <- cpp_list_assign
#' @export
#' @rdname lists
list_modify <- list_assign

#' @export
#' @rdname lists
list_combine <- function(..., .args = NULL){
  cpp_list_c(.Call(`_cheapr_cpp_list_args`, list(...), .args))
}

#' @export
#' @rdname lists
list_drop_null <- cpp_drop_null

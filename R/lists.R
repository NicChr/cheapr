#' List utilities
#'
#' @description
#' Functions to help work with lists.
#'
#' @param x A list.
#' @param length Length of list.
#' @param default Default value for each list element.
#'
#' @returns
#' `lengths_()` returns the list lengths. \cr
#' `unlisted_length()` is an alternative to `length(unlist(x))`. \cr
#' `new_list()` is like `vector("list", length)` but also allows you to specify
#' a default value for each list element. This can be useful for
#' initialising with a catch-all value so that when you unlist you're guaranteed
#' a list of length >= to the specified length.
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
#' unlisted_length(l) - num_na(l) # Number of non-NA elements
#'
#' # We can create and initialise a new list with a default value
#' l <- new_list(20, 0L)
#' l[1:5]
#' # This works well with vctrs_list_of objects
#' vctrs::new_list_of(l, ptype = integer())[1:5]
#' @export
#' @rdname lists
lengths_ <- cpp_lengths
#' @export
#' @rdname lists
unlisted_length <- cpp_r_unnested_length
#' @export
#' @rdname lists
new_list <- function(length = 0L, default = NULL){
  if (base::length(length) != 1){
    stop("length must be a vector of length 1")
  }
  cpp_new_list(as.numeric(length), default)
}


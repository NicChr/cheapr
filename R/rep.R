#' cheapr style repeat functions
#'
#' @name rep
#'
#' @param x A vector or data frame.
#' @param times `[integer(n)]` A vector of times to repeat elements of `x`. Can be length 1
#' or the same length as `vector_length(x)`.
#' @param length `[integer(1)]` - Length of the recycled result.
#' @param each `[integer(1)]` - How many times to repeat out each element of `x`.
#'
#' @returns
#' Repeated out object.
#'
#' @rdname rep
#' @export
cheapr_rep <- cpp_rep
#' @rdname rep
#' @export
rep_ <- cheapr_rep
#' @rdname rep
#' @export
cheapr_rep_len <- cpp_rep_len
#' @rdname rep
#' @export
rep_len_ <- cheapr_rep_len
#' @rdname rep
#' @export
cheapr_rep_each <- cpp_rep_each
#' @rdname rep
#' @export
rep_each_ <- cheapr_rep_each

#' Recycle objects to a common size
#'
#' @description
#' A convenience function to recycle R objects to either a common or specified
#' size.
#'
#' @param ... Objects to recycle.
#' @param length Optional length to recycle objects to.
#' @param .args An alternative to `...` so you can supply arguments directly
#' in a list. \cr
#' This is equivalent to `do.call(f, .args)` but much more efficient.
#'
#' @returns
#' A list of recycled R objects.
#'
#' @details
#' Data frames are recycled by recycling their rows. \cr
#' `recycle()` is optimised to only recycle objects that need recycling. \cr
#' `NULL` objects are ignored and not recycled or returned.
#'
#' @examples
#' library(cheapr)
#'
#' # Recycles both to size 10
#' recycle(Sys.Date(), 1:10)
#'
#' # Any vectors of zero-length are all recycled to zero-length
#' recycle(integer(), 1:10)
#'
#' # Unless length is supplied
#' recycle(integer(), 1:10, length = 10)
#'
#' # Data frame rows are recycled
#' recycle(sset(iris, 1:3), length = 9)
#'
#' # To recycle objects in a list, use `.args`
#' my_list <- list(from = 1L, to = 10L, by = seq(0.1, 1, 0.1))
#' recycle(.args = my_list)
#'
#' @export
recycle <- function (..., length = NULL, .args = NULL){
  .Call(`_cheapr_cpp_recycle`, .Call(`_cheapr_cpp_list_args`, list(...), .args), length)
}

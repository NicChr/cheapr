#' Recycle objects to a common size
#'
#' @description
#' A convenience function to recycle R objects to either a common or specified
#' size.
#'
#' @param ... Objects to recycle.
#' @param length Optional length to recycle objects to.
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
#' recycle(Sys.Date(), 1:10)
#'
#' # Any vectors of zero-length are all recycled to zero-length
#' recycle(integer(), 1:10)
#'
#' # Data frame rows are recycled
#' recycle(sset(iris, 1:3), length = 3 * 3)
#'
#' # To recycle list items, use `do.call()`
#' my_list <- list(from = 1L, to = 10L, by = seq(0.1, 1, 0.1))
#' do.call(recycle, my_list)
#'
#' @export
recycle <- function (..., length = NULL){
  out <- cpp_list_rm_null(list(...))
  lens <- lengths_(out)
  if (is.null(length)) {
    if (length(lens)) {
      N <- max(lens)
    }
    else {
      N <- 0L
    }
  }
  else {
    N <- length
  }
  N <- N * (!collapse::anyv(lens, 0L))
  recyclei <- which_(lens != N)
  if (length(recyclei)){
    out[recyclei] <- lapply(out[recyclei], cheapr_rep_len, N)
  }
  out
}

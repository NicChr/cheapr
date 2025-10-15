#' A cheapr version of `c()`
#'
#' cheapr's version of `c()`. It is quite a bit faster for atomic vectors
#' and combines data frame rows instead of cols.
#'
#' @param ... Objects to combine.
#' @param .args An alternative to `...` so you can supply arguments directly
#' in a list. \cr
#' This is equivalent to `do.call(f, .args)` but much more efficient.
#'
#' @returns
#' Combined objects.
#'
#' @examples
#' library(cheapr)
#'
#' # Combine just like `c()`
#' c_(1, 2, 3:5)
#'
#' # It combines rows by default instead of cols
#' c_(new_df(x = 1:3), new_df(x = 4:10))
#'
#' # If you have a list of objects you want to combine
#' # use `.args` instead of `do.call` as it's more efficient
#'
#' list_of_objs <- rep_(list(0), 10^4)
#'
#' bench::mark(
#'   do.call(c, list_of_objs),
#'   do.call(c_, list_of_objs),
#'   c_(.args = list_of_objs) # Fastest
#' )
#'
#' @rdname cheapr_c
#' @export
cheapr_c <- function(..., .args = NULL){
  .Call(`_cheapr_cpp_c`, .Call(`_cheapr_cpp_list_args`, list(...), .args))
}
#' @rdname cheapr_c
#' @export
c_ <- cheapr_c

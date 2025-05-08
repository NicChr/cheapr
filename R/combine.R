#' A cheapr version of `c()`
#'
#' cheapr's version of `c()`. It is quite a bit faster for atomic vectors
#' and combines data frame rows instead of cols.
#'
#' @param ... Objects to combine.
#' @param .args An alternative to `...` for easier programming with lists.
#'
#' @returns
#' Combined objects.
#'
#' @examples
#' library(cheapr)
#'
#' # Combine just like `c()`
#' cheapr_c(1, 2, 3:5)
#'
#' # It combines rows by default instead of cols
#' cheapr_c(new_df(x = 1:3), new_df(x = 4:10))
#'
#' # If you have a list of objects you want to combine
#' # use `.args` instead of `do.call` as it's more efficient
#'
#' list_of_objs <- rep(list(0), 10^4)
#'
#'  bench::mark(
#'     do.call(cheapr_c, list_of_objs),
#'     cheapr_c(.args = list_of_objs)
#'   )
#' @export
cheapr_c <- function(..., .args = NULL){
  .Call(`_cheapr_cpp_c`, .Call(`_cheapr_cpp_list_args`, list(...), .args))
}

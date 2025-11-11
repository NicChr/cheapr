#' Switch between dot-dot-dot and a list of args
#'
#' @description
#' `switch_args()` is primarily used as a helper when writing functions that
#' use the dot-dot-dot argument `...`. \cr
#' cheapr uses it frequently to give more freedom to the user in
#' how they pass arguments to
#' functions that take a variable number of arguments.\cr
#' See examples for info.
#'
#' @param ... Arguments passed individually.
#' @param .args Alternative list of arguments.
#' Either `...` or `.args` can be used, not both.
#'
#' @details
#' Using `switch_args` simply allows the user to avoid having to use
#' `do.call(fn, args)`. This can be advantageous for developers who write
#' compiled (C/C++) functions that accept lists of arguments.
#' cheapr internally uses this framework for performance critical
#' functions such as `cheapr::c_()` which internally calls
#' `cheapr:::cpp_c()`, a function that takes one list of vectors and
#' combines them into one vector. The equivalent of `cheapr::c_(.args = args)`
#' would be the less efficient `do.call(cheapr::c_, args)`.
#'
#' @returns
#' A list of arguments
#'
#' @examples
#' library(cheapr)
#'
#' # Flexibly create a data frame
#' base_new_df <- function(..., .args = NULL){
#'
#'   args <- switch_args(..., .args = .args)
#'
#'   list2DF(args)
#' }
#'
#' # Normal usage
#' base_new_df(x = 1, y = 2)
#'
#' # Alternatively supplying a list of args instead
#' base_new_df(.args = list(x = 1, y = 2))
#'
#' # cheapr::new_df does something similar
#' new_df(x = 1, y = 2)
#' new_df(.args = list(x = 1, y = 2))
#'
#' @export
switch_args <- function(..., .args = NULL){
  .Call(`_cheapr_cpp_list_args`, list(...), .args)
}

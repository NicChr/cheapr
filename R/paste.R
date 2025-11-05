#' Fast string concatenation using C++
#'
#' @rdname paste
#' @param ... Character vectors to concatenate.
#'
#' @param .sep `[character(1)]` - A string to separate the supplied strings.
#' @param .args An alternative to `...` so you can supply arguments directly
#' in a list. \cr
#' This is equivalent to `do.call(f, .args)` but much more efficient.
#'
#' @examples
#' library(cheapr)
#'
#' # Normal usage
#' paste_("Hello", "and", "Goodbye", sep = " ")
#' paste0_(100, "%")
#'
#' # Recycling with zero-length vectors
#' paste_("hello", character(), letters)
#'
#' # Concatenating many objects is very fast via `.args`
#' strings <- sample_(letters, 5e04, TRUE) |>
#'   with_local_seed(1) |>
#'   as.list()
#' mark(paste_(.args = strings), do.call(paste, strings))
#'
#' @export
paste_ <- function(..., .sep = " ", .args = NULL){
  .Call(`_cheapr_cpp_paste`, .Call(`_cheapr_cpp_list_args`, list(...), .args), .sep)
}
#' @rdname paste
#' @export
paste0_ <- function(..., .args = NULL){
  .Call(`_cheapr_cpp_paste`, .Call(`_cheapr_cpp_list_args`, list(...), .args), "")
}

#' Fast string concatenation using C++
#'
#' @name strings
#'
#' @param ... Character vectors to concatenate.
#'
#' @param sep `[character(1)]` - A string to separate the supplied strings.
#' @param collapse Optional string to collapse concatenated strings into
#' one string (character vector of length 1).
#' @param .args An alternative to `...` so you can supply arguments directly
#' in a list. \cr
#' This is equivalent to `do.call(f, .args)` but much more efficient.
#'
#' @examples
#' library(cheapr)
#' # Normal usage
#' paste_("Hello", "and", "Goodbye", sep = " ")
#' paste_(100, "%")
#'
#' paste_(letters, LETTERS)
#' paste_(letters, LETTERS, collapse = "")
#'
#' # Both concatenating and collapsing
#' paste_(letters, LETTERS, sep = ",", collapse = " next letter ")
#' # This is the same as above
#' paste_(letters, LETTERS, collapse = ", next letter ")
#'
#' # Recycling with zero-length vectors
#' paste_("hello", character(), letters)
#' paste_("hello", character(), letters, collapse = "")
#'
#' # Concatenating many objects is very fast via `.args`
#'
#' sampled_letters <- sample_(letters, 5e04, TRUE)
#'
#' # Pasting multiple character vectors
#' mark(
#'   paste_(sampled_letters, sep = ","),
#'   paste(sampled_letters, sep = ",")
#' )
#' # Collapsing is very fast compared to base R
#' mark(
#'   paste_(sampled_letters, collapse = ""),
#'   paste(sampled_letters, collapse = "")
#' )
#'
#' strings <- sampled_letters |>
#'   with_local_seed(1) |>
#'   as.list()
#'
#' strings <- lapply(strings, rep_len_, 3)
#'
#' library(bench)
#' mark(
#'   paste_(.args = strings),
#'   do.call(paste0, strings)
#' )
#' mark(
#'   paste_(.args = strings, collapse = ""),
#'   do.call(paste_, c(strings, list(collapse = ""))),
#'   do.call(paste0, c(strings, list(collapse = "")))
#' )
#'
#' @rdname strings
#' @export
paste_ <- function(..., sep = "", collapse = NULL, .args = NULL){
  .Call(
    `_cheapr_cpp_paste`, .Call(`_cheapr_cpp_list_args`, list(...), .args),
    sep, collapse
  )
}

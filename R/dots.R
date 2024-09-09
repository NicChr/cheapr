#' Turn dot-dot-dot (`...`) into a named list
#'
#' @description
#' A fast and useful function for always returning a named list from `...`
#'
#'
#' @param ... Key-value pairs.
#'
#' @returns A named list.
#'
#' @rdname dots
#' @export
named_dots <- function(...){
  dots <- list(...)

  dot_nms <- names(dots)

  if (is.null(dot_nms)){
    names(dots) <- expr_names(...)
  } else if (!all(nzchar(dot_nms))){
    empty <- which_(!nzchar(dot_nms))
    expr_nms <- expr_names(...)
    dot_nms[empty] <- expr_nms[empty]
    names(dots) <- dot_nms
  }
  dots
}
expr_names <- function(...){
  as.character(substitute(c(...))[-1L])
  # vapply(substitute(alist(...))[-1L], deparse2, "", USE.NAMES = FALSE)
}
list_named <- function(...){
  cpp_list_rm_null(named_dots(...))
}

#' Turn dot-dot-dot (`...`) into a named list
#'
#' @description
#' A fast and useful function for always returning a named list from `...`
#'
#'
#' @param ... Key-value pairs.
#' @param .keep_null Should `NULL` entries be kept? Default is `TRUE`.
#'
#' @returns A named list.
#'
#' @rdname dots
#' @export
named_list <- function(..., .keep_null = TRUE){
  dots <- list(...)
  dot_nms <- names(dots)

  if (is.null(dot_nms)){
    names(dots) <- expr_names(...)
  } else if (!all(nzchar(dot_nms))){
    empty <- which_(nzchar(dot_nms), invert = TRUE)
    expr_nms <- expr_names(...)
    dot_nms[empty] <- expr_nms[empty]
    names(dots) <- dot_nms
  }
  if (!.keep_null){
    dots <- cpp_list_rm_null(dots)
  }
  dots
}
expr_names <- function(...){
  # as.character(substitute(c(...))[-1L])
  vapply(exprs(...), deparse2, "", USE.NAMES = FALSE)
}

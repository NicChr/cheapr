#' Fast casting/coercing of R objects
#'
#' @description
#' `cast()` is type-commutative, meaning the order `x` and `y` doesn't
#' affect the outcome type.
#' `cast()` will attempt to cast `x` into a common type between `x` and `y`.
#'
#' @param x A vector.
#' @param y A vector.
#' @param ... Vectors.
#' @param .args An alternative to `...` so you can supply arguments directly
#' in a list. \cr
#' This is equivalent to `do.call(f, .args)` but much more efficient.
#'
#' @returns
#' `cast()` coerces `x` into a common type between `x` and `y`.
#' `cast_common()` coerces all supplied vectors into a common type between them.
#'
#' @rdname cast
#' @export
cast <- function(x, y) {
  .Call(`_cheapr_cpp_cast`, x, y)
}

#' @rdname cast
#' @export
cast_common <- function(..., .args = NULL){
  .Call(`_cheapr_cpp_cast_all`, .Call(`_cheapr_cpp_list_args`, list(...), .args))
}

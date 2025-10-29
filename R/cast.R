#' Fast casting/coercing of R objects
#'
#' @description
#' `cast_common()` is type-commutative, meaning the order of objects doesn't
#' affect the outcome type.
#' `cast()` will attempt to cast `x` into an object similar to `y`.
#'
#' @param x A vector.
#' @param y A vector.
#' @param ... Vectors.
#' @param .args An alternative to `...` so you can supply arguments directly
#' in a list. \cr
#' This is equivalent to `do.call(f, .args)` but much more efficient.
#'
#' @returns
#' `cast()` will attempt to cast `x` into an object similar to `y`.
#' `cast_common()` coerces all supplied vectors into a common type between them.
#'
#' @export
cast <- function(x, y){
  .Call(`_cheapr_cpp_cast`, x, y)
}

#' @rdname cast
#' @export
cast_common <- function(..., .args = NULL){
  .Call(`_cheapr_cpp_cast_common`, .Call(`_cheapr_cpp_list_args`, list(...), .args))
}

base_cast <- function(x, template){
  if (is.null(template)){
    x
  } else {
    x[0] <- template[0]
    x
    # c(x, template[0])
  }
}

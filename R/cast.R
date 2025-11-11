#' Fast casting/coercing of R objects
#'
#' @description
#' `cast_common()` is type-commutative, meaning the order of objects doesn't
#' affect the outcome type.
#' `cast()` will attempt to cast `x` into an object similar to `archetype`.
#'
#' @param x A vector.
#' @param archetype An archetype vector.
#' @param ... Vectors.
#' @param .args An alternative to `...` so you can supply arguments directly
#' in a list. \cr
#' This is equivalent to `do.call(f, .args)` but much more efficient.
#'
#' @returns
#' `cast()` will attempt to cast `x` into an object similar to `archetype`. \cr
#' `cast_common()` coerces all supplied vectors into a
#' common type between them. \cr
#' `archetype()` returns the zero-length template/archetype of `x`. \cr
#' `archetype_common()` returns the common zero-length
#' template between all supplied vectors. \cr
#' `r_type()` will return the internal cheapr-defined type of `x` as a
#' character vector of length 1. This will usually match `class(x)`
#' but not always. \cr
#' `r_type_common()` returns the common type between all objects.
#'
#' @rdname cast
#' @export
cast <- function(x, archetype){
  .Call(`_cheapr_cpp_cast`, x, archetype)
}

#' @rdname cast
#' @export
cast_common <- function(..., .args = NULL){
  .Call(`_cheapr_cpp_cast_common`, .Call(`_cheapr_cpp_list_args`, list(...), .args))
}

#' @rdname cast
#' @export
archetype <- function(x){
  .Call(`_cheapr_cpp_common_template`, .args = list(x))
}

#' @rdname cast
#' @export
archetype_common <- function(..., .args = NULL){
  .Call(`_cheapr_cpp_common_template`, .Call(`_cheapr_cpp_list_args`, list(...), .args))
}

#' @rdname cast
#' @export
r_type <- function(x) {
  .Call(`_cheapr_cpp_type`, x)
}

#' @rdname cast
#' @export
r_type_common <- function(..., .args = NULL){
  .Call(
    `_cheapr_cpp_type`,
    .Call(`_cheapr_cpp_common_template`, .Call(`_cheapr_cpp_list_args`, list(...), .args))
  )
}

# Internal function, do not use
base_cast <- function(x, template){
  if (is.null(template)){
    NULL
  } else if (is.null(x)){
    template[0]
  } else {
    x[0] <- template[0]
    x
    # c(template[0], x)
  }
}

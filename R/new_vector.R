#' Fast multi-threaded vector initialisation
#'
#' @param n `[integer(1)]` - Length of vector.
#' @param names `[character(n)]` - Names of initialised vector.
#' @param default Default value to initialise the vector with.
#'
#' @returns
#' New vector.
#'
#' @seealso [new_list]
#'
#' @rdname new_vector
#' @export
new_logical <- function(n = 0L, names = NULL, default = FALSE){
  if (length(default) != 1){
    stop("`default` must be a length one vector")
  }
  attrs_modify(
    rep_len_(as.logical(default), n),
    names = names,
    .set = TRUE
  )
}
#' @rdname new_vector
#' @export
new_integer <- function(n = 0L, names = NULL, default = 0L){
  if (length(default) != 1){
    stop("`default` must be a length one vector")
  }
  attrs_modify(
    rep_len_(as.integer(default), n),
    names = names,
    .set = TRUE
  )
}
#' @rdname new_vector
#' @export
new_double <- function(n = 0L, names = NULL, default = 0){
  if (length(default) != 1){
    stop("`default` must be a length one vector")
  }
  attrs_modify(
    rep_len_(as.double(default), n),
    names = names,
    .set = TRUE
  )
}
#' @rdname new_vector
#' @export
new_character <- function(n = 0L, names = NULL, default = ""){
  if (length(default) != 1){
    stop("`default` must be a length one vector")
  }
  attrs_modify(
    rep_len_(as.character(default), n),
    names = names,
    .set = TRUE
  )
}
#' @rdname new_vector
#' @export
new_complex <- function(n = 0L, names = NULL, default = 0i){
  attrs_modify(
    rep_len_(as.complex(default), n),
    names = names,
    .set = TRUE
  )
}
#' @rdname new_vector
#' @export
new_raw <- function(n = 0L, names = NULL, default = as.raw(0x00)){
  if (length(default) != 1){
    stop("`default` must be a length one vector")
  }
  attrs_modify(
    rep_len_(as.raw(default), n),
    names = names,
    .set = TRUE
  )
}

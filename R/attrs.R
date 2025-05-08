#' Add and remove attributes
#'
#' @description
#' Simple tools to add and remove attributes, both normally and in-place.
#' To remove specific attributes, set those attributes to `NULL`.
#'
#'
#' @param x Object to add/remove attributes.
#' @param ... Named attributes, e.g 'key = value'.
#' @param .set Should attributes be added in-place without shallow-copying `x`?
#' Default is `FALSE`.
#' @param .args An alternative to `...` for easier programming with lists.
#'
#' @seealso [shallow_copy]
#'
#' @name attrs
#'
#' @returns
#' The object `x` with attributes removed or added.
#'
#' @rdname attrs
#' @export
attrs_add <- function(x, ..., .set = FALSE, .args = NULL){
  .Call(
    `_cheapr_cpp_set_add_attributes`,
    if (.set) x else .Call(`_cheapr_cpp_shallow_copy`, x),
    .Call(`_cheapr_cpp_list_args`, list(...), .args), TRUE
  )
}
#' @rdname attrs
#' @export
attrs_rm <- function(x, .set = FALSE){
  .Call(`_cheapr_cpp_set_rm_attributes`, if (.set) x else .Call(`_cheapr_cpp_shallow_copy`, x))
}

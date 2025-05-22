#' Rebuild an object from a template
#'
#' @param x An object in which carefully selected attributes
#' will be copied into from `template`.
#' @param template A template object used to copy attributes into `x`.
#' @param shallow_copy Should `x` be shallow copied before rebuilding?
#' Default is `TRUE`.
#' @param ... Further arguments passed onto methods.
#'
#' @details
#' In R attributes are difficult to work with. One big reason for this is
#' that attributes may or may not be independent of the data.
#' Date vectors for example have attributes completely independent of the data
#' and hence if the attributes are removed at any point, they can easily be
#' re-added without any calculations. Factors have almost data-independent
#' attributes with an exception being when factors are combined.
#' In some cases it is not possible to rebuild attributes from the data
#' alone.
#'
#' You can add your own `rebuild` method for an object not covered
#' by the methods here.
#'
#' @returns
#' An object similar to `template`.
#'
#' @rdname rebuild
#' @export
rebuild <- function(x, template, ...){
  UseMethod("rebuild", template)
}
# rebuild.default <- function(x, template, shallow_copy = TRUE, ...){
#   cpp_reconstruct(
#     x, template, c("names", "dim", "dimnames", "row.names", "tsp", "comment"),
#     cpp_setdiff(
#       names(attributes(template)),
#       c("names", "dim", "dimnames", "row.names", "tsp", "comment")
#     ), shallow_copy
#   )
# }
#' @rdname rebuild
#' @export
rebuild.data.frame <- function(x, template, shallow_copy = TRUE, ...){
  .Call(`_cheapr_cpp_rebuild`, x, template, c("names", "row.names"), "class", shallow_copy)
}
#' @rdname rebuild
#' @export
rebuild.data.table <- function(x, template, shallow_copy = TRUE, ...){
  collapse::qDT(
    cpp_rebuild(x, template, c("names", "row.names", "sorted"), "class", shallow_copy),
    class = class(template)
  )
}


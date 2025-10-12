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
#' @rdname rebuild
#' @export
rebuild.data.frame <- function(x, template, shallow_copy = TRUE, ...){

  names <- names(x)
  row_names <- .set_row_names(length(attr(x, "row.names", TRUE)))

  attrs_modify(
    attrs_clear(x, .set = !shallow_copy),
    names = names, row.names = row_names,
    class = "data.frame", .set = TRUE
  )
}
#' @rdname rebuild
#' @export
rebuild.data.table <- function(x, template, shallow_copy = TRUE, ...){

  names <- names(x)
  row_names <- .set_row_names(length(attr(x, "row.names", TRUE)))
  sorted <- attr(x, "sorted", TRUE)

  # qDT() will internally add a true length
  collapse::qDT(
    attrs_modify(
      attrs_clear(x, .set = !shallow_copy),
      names = names, row.names = row_names, sorted = sorted,
      class = c("data.table", "data.frame"),
      .set = TRUE
    )
  )
}
#' @rdname rebuild
#' @export
rebuild.tbl_df <- function(x, template, shallow_copy = TRUE, ...){

  names <- names(x)
  row_names <- .set_row_names(length(attr(x, "row.names", TRUE)))

  attrs_modify(
    attrs_clear(x, .set = !shallow_copy),
    names = names, row.names = row_names,
    class = c("tbl_df", "tbl", "data.frame"),
    .set = TRUE
  )
}

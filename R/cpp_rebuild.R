#' Low-level attribute re-constructor
#'
#' @param target Target object you wish to rebuild attributes on.
#' @param source Source object to copy attributes from.
#' @param target_attr_names `[character(n)]` Names of target attributes to keep.
#' @param source_attr_names `[character(n)]` Names of source attributes
#' to copy onto target.
#' @param shallow_copy `[logical(1)]` Should target be shallow copied before re-building?
#' If `FALSE` attributes are added in-place.
#'
#' @details
#' `cpp_rebuild()` is mostly a convenience function to help with choosing
#' exactly which attributes to copy onto the target object.
#' `rebuild()` is a related generic function with rebuild methods for
#' common objects (currently only `tbl_df`, `data.frame` and `data.table`).
#' For examples of further rebuild methods, see the fastplyr package.
#'
#' To modify attributes yourself you can of course use base R attribute functions
#' like `attr()` and `attributes()` or cheapr's more convenient `attrs_modify`.
#'
#' @seealso [rebuild] [attrs_modify]
#'
#' @returns
#' An object similar to `source`.
#'
#' @export
cpp_rebuild <- cpp_rebuild

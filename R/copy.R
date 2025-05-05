#' Copy R objects
#'
#' @description
#' `shallow_copy()` and `deep_copy()` are just wrappers to the
#' R C API functions `Rf_shallow_duplicate()` and `Rf_duplicate()`
#' respectively. `semi_copy()` is something in between whereby it
#' fully copies the data but only shallow copies the attributes.
#'
#' @details
#' Shallow duplicates are mainly useful for adding attributes to objects
#' in-place as well assigning vectors to shallow copied lists in-place.
#'
#' Deep copies are generally useful for ensuring an object is fully
#' duplicated, including all attributes associated with it.
#' Deep copies are generally expensive and should be
#' used with care.
#'
#' To summarise:
#'
#' * `shallow_copy` - Shallow copies data and attributes
#' * `semi_copy` - Deep copies data and shallow copies attributes
#' * `deep_copy` - Deep copies both data and attributes
#'
#' It is recommended to use these functions only if you know what you are doing.
#'
#' @param x An object to shallow, semi, or deep copy.
#'
#' @returns
#' A shallow, semi or deep copied R object.
#'
#' @name copy
#'
#' @examples
#'
#' library(cheapr)
#' library(bench)
#' df <- new_df(x = sample.int(10^4))
#'
#' # Note the memory allocation
#' mark(shallow_copy(df), iterations = 1)
#' mark(deep_copy(df), iterations = 1)
#'
#' # In both cases the address of df changes
#'
#' address(df);address(shallow_copy(df));address(deep_copy(df))
#'
#' # When shallow-copying attributes are not duplicated
#'
#' address(attr(df, "names"));address(attr(shallow_copy(df), "names"))
#'
#' # They are when deep-copying
#'
#' address(attr(df, "names"));address(attr(deep_copy(df), "names"))
#'
#' # Adding an attribute in place with and without shallow copy
#' invisible(attrs_add(df, key = TRUE, .set = TRUE))
#' attr(df, "key")
#'
#' # Remove attribute in-place
#' invisible(attrs_add(df, key = NULL, .set = TRUE))
#'
#' # With shallow copy
#' invisible(attrs_add(shallow_copy(df), key = TRUE, .set = TRUE))
#'
#' # 'key' attr was only added to the shallow copy, and not the original df
#' attr(df, "key")
#'
#' @rdname copy
#' @export
shallow_copy <- function(x){
  .Call(`_cheapr_cpp_shallow_copy`, x)
}
#' @rdname copy
#' @export
semi_copy <- function(x){
  .Call(`_cheapr_cpp_semi_copy`, x)
}
#' @rdname copy
#' @export
deep_copy <- function(x){
  .Call(`_cheapr_r_copy`, x)
}

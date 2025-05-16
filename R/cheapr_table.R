#' Fast frequency tables - Still experimental
#'
#' @description
#' This is not a one-to-one copy of `base::table()` as some behaviours differ.
#' It is more flexible as it accepts inputs such as data frames and
#' `vctrs_rcrd` objects.
#'
#'
#' @param ... `>=1` objects that can be converted to a factor through
#' `cheapr::factor_()`.
#' @param names Should level names be kept? Default is `TRUE`.
#' @param order Should result be ordered by level names? Default is
#' `FALSE`.
#' @param na_exclude Should `NA` values be excluded? Default is `FALSE`.
#' @param classed Should a `table` object be returned? Default is `FALSE`
#'
#' @returns
#' A named integer vector if one object is supplied, otherwise an
#' array.
#'
#' @export
cheapr_table <- function(..., names = TRUE, order = FALSE, na_exclude = FALSE,
                         classed = FALSE){

  vecs <- list(...)

  to_factor <- function(x, na.excl){
    if (is.factor(x)){
      if (na_exclude){
        levels_drop_na(x)
      } else if (any_na(x)){
        levels_add_na(x)
      } else {
        x
      }
    } else {
      factor_(x, order = order, na_exclude = na_exclude)
    }
  }

  factors <- lapply(vecs, to_factor, na_exclude)

  if (length(factors) == 1){
    f <- factors[[1L]]
    lvls <- attr(f, "levels")
    out <- tabulate(f, nbins = length(lvls))
    if (names){
      dim_names <- list(lvls)
      names(dim_names) <- names(vecs)
      names(out) <- lvls
    } else {
      dim_names <- NULL
    }
    if (classed){
      out <- array(out, length(out), dim_names)
      class(out) <- "table"
    }
  } else {
    out <- do.call(collapse::qtab, c(factors, list(dnn = NULL)), envir = parent.frame())
    if (!classed){
      attrs_add(out, class = NULL, sorted = NULL, weighted = NULL, .set = TRUE)
    }
  }
  out
}

# counts <- function(x, sort = FALSE){
#   if (!cpp_is_simple_atomic_vec(x)){
#     stop("`x` must be an atomic vector")
#   }
#   if (sort){
#     groups <- collapse::GRP(
#       x, sort = TRUE, return.order = FALSE, return.groups = TRUE
#     )
#     n_groups <- groups[["N.groups"]]
#     starts <- groups[["group.starts"]]
#     ids <- groups[["group.id"]]
#     keys <- groups[["groups"]][[1L]]
#   } else {
#     groups <- collapse::group(x, starts = TRUE)
#     ids <- groups
#     n_groups <- attr(groups, "N.groups", TRUE)
#     starts <- attr(groups, "starts", TRUE)
#     keys <- cheapr::sset(x, starts)
#   }
#   new_df(
#     key = keys,
#     count = cpp_tabulate(ids, n_groups)
#   )
# }

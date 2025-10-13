#' An alternative `unique` function
#'
#' @description
#' `unique_()` is a usually faster alternative to `unique()` with optional
#' sorting included. The internal API of this function is designed to be
#' simple and generic to allow for working with all kinds of objects that can
#' be reduced to a unique set.
#'
#' Internally `unique_()` calculates unique group IDs for the given vector
#' in the range `[1, n]` where `1` denotes the first group and `n` denotes
#' the nth group.
#' This function will work correctly as long as there
#' is a correctly implemented `collapse::GRP` method and a `[` method
#' for the object.
#' In the future cheapr will include a `group_id` S3 generic to replace the use
#' of `collapse::GRP` here, of which is arguably more difficult to write correct
#' methods for.
#'
#' @param x A vector (or data frame).
#' @param sort Should unique result be sorted? Default is `FALSE`.
#'
#' @returns
#' A unique vector (or data frame).
#'
#' @examples
#'
#' library(cheapr)
#'
#' x <- rep_(3:1, 3)
#' unique_(x)
#' unique_(x, sort = TRUE)
#'
#' # Unique rows
#' iris |>
#'   sset(j = c("Petal.Width", "Species")) |>
#'   unique_()
#' @export
unique_ <- function(x, sort = FALSE){

  # In the future, `group_id()` will be moved from fastplyr
  # to cheapr to make things simpler
  groups <- collapse::GRP(
    x, sort = sort,
    return.groups = FALSE,
    return.order = sort
  )
  group_ids <- groups[["group.id"]]
  n_groups <- groups[["N.groups"]]

  # If `x` is already unique, return x
  if (!sort && vector_length(x) == n_groups){
    return(x)
  }
  if (sort && isTRUE(attr(groups[["order"]], "sorted"))){
    return(x)
  }

  start_locs <- cpp_group_starts(group_ids, n_groups)
  sset(x, start_locs)
}

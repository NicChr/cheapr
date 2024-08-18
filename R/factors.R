#' A faster version of `factor()`
#'
#' @description
#' A fast version of `factor()` using the collapse package. \cr
#'
#' There are some additional utilities such as
#' `levels_factor()` which returns the levels of a factor, as a factor,
#' `used_levels()` which returns the used levels of a factor,
#' `unused_levels()` which returns the unused levels of a factor,
#' `add_na_level()` which adds an explicit `NA` level,
#' `drop_na_level()` which drops the `NA` level,
#' `drop_levels()` which drops unused factor levels,
#' and finally `reorder_levels()` which reorders the levels of `x`
#' based on `y` using the ordered median values of `y` for each level.
#'
#' @returns
#' A `factor` or `character` in the case of `used_levels` and `unused_levels`.
#'
#' @param x A vector.
#' @param levels Optional factor levels.
#' @param order Should factor levels be sorted? Default is `TRUE`.
#' It typically is faster to set this to `FALSE`, in which case the levels
#' are sorted by order of first appearance.
#' @param na_exclude Should `NA` values be excluded from the factor levels?
#' Default is `TRUE`.
#' @param ordered Should the result be an ordered factor?
#' @param name Name of `NA` level.
#' @param where Where should `NA` level be placed? Either first or last.
#' @param order_by A vector to order the levels of `x` by using the medians of
#' `order_by`.
#' @param decreasing Should the reordered levels be in decreasing order?
#' Default is `FALSE`.
#'
#' @details
#' This operates similarly to `collapse::qF()`. \cr
#' The main difference internally is that `collapse::funique()` is used
#' and therefore s3 methods can be written for it. \cr
#' Furthermore, for date-times `factor_` differs in that it differentiates
#' all instances in time whereas `factor` differentiates calendar times.
#' Using a daylight savings example where the clocks go back: \cr
#' `factor(as.POSIXct(1729984360, tz = "Europe/London") + 3600 *(1:5))`
#' produces 4 levels whereas \cr
#' `factor_(as.POSIXct(1729984360, tz = "Europe/London") + 3600 *(1:5))`
#' produces 5 levels.
#'
#' @export
#' @rdname factors
factor_ <- function(
    x = integer(), levels = NULL, order = TRUE,
    na_exclude = TRUE,
    ordered = is.ordered(x)
){
  if (is.null(x)){
    x <- integer()
  }
  if (is.null(levels)){
    lvls <- collapse::funique(x, sort = order, na.last = TRUE)
  } else {
    lvls <- levels
  }
  if (na_exclude && any_na(lvls)){
    if (order && is.null(levels)){
      lvls <- sset(lvls, seq_len(cpp_vec_length(lvls) - 1L))
    } else {
      lvls <- na_rm(lvls)
    }
  }
  out <- collapse::fmatch(x, lvls, overid = 2L)
  if (inherits(lvls, "data.frame")){
    fct_lvls <- do.call(paste, c(lvls, list(sep = "_")))
  } else {
    fct_lvls <- as.character(lvls)
  }
  if (inherits(x, "POSIXt") && collapse::any_duplicated(fct_lvls)){
    fct_lvls <- paste(fct_lvls, as.POSIXlt(lvls)$zone)
  }
  # if (!identical(levels, labels)){
  #   fct_lvls <- as.character(labels)
  # }
  attr(out, "levels") <- fct_lvls
  class(out) <- c(if (ordered) "ordered" else character(), "factor")
  # class(out) <- c("cheapr_factor", if (ordered) "ordered" else character(),
  #                 # if (!na_exclude) "na_included" else character(),
  #                 "factor")
  # attr(out, "data") <- x
  out
}
# print.cheapr_factor <- function(x, ...){
#   class(x) <- setdiff(class(x), "cheapr_factor")
#   attr(x, "data") <- NULL
#   # class(x) <- setdiff(class(x), "cheapr_factor")
#   # NextMethod("print")
#   print(x, ...)
# }
# as.character.cheapr_factor <- function(x, ...){
#   data <- attr(x, "data")
#   if (!is.null(data)){
#     x <- data
#   }
#   as.character(x, ...)
#   # x <- attr(x, "data")
#   # class(x) <- setdiff(class(x), "cheapr_factor")
#   # NextMethod("as.character")
# }
# `[.cheapr_factor` <- function(x, ...){
#   data <- attr(x, "data")
#   attr(x, "data") <- NULL
#   out <- NextMethod("[")
#   if (!is.null(data)){
#     attr(out, "data") <- data[...]
#   }
#   out
# }
#
# as.double.cheapr_factor <- function(x, ...){
#   data <- attr(x, "data")
#   if (!is.null(data)){
#     x <- data
#   }
#   as.double(x, ...)
#   # NextMethod("as.numeric")
# }
# as.integer.cheapr_factor <- function(x, ...){
#   data <- attr(x, "data")
#   if (!is.null(data)){
#     x <- data
#   }
#   as.integer(x, ...)
#   # NextMethod("as.integer")
# }
#' @export
#' @rdname factors
levels_factor <- function(x){
  check_is_factor(x)
  lvls <- levels(x)
  out <- seq_along(lvls)
  attr(out, "levels") <- lvls
  class(out) <- class(x)
  out
}
check_is_factor <- function(x){
  if (!is.factor(x)){
    stop("x must be a factor")
  }
}
factor_as_character <- function(x){
  levels(x)[unclass(x)]
}
which_used_levels <- function(x){
  which_in(levels_factor(x), x)
}
which_unused_levels <- function(x){
  which_not_in(levels_factor(x), x)
}
#' @export
#' @rdname factors
used_levels <- function(x){
  levels(x)[which_used_levels(x)]
  # levels(x)[which_in(seq_along(levels(x)), unclass(x))]
  # factor_as_character(intersect_(levels_factor(x), x))
}
#' @export
#' @rdname factors
unused_levels <- function(x){
  levels(x)[which_unused_levels(x)]
  # levels(x)[which_not_in(seq_along(levels(x)), unclass(x))]
  # factor_as_character(setdiff_(levels_factor(x), x))
  # as.character(setdiff_(levels_factor(x), x))
}
#' @export
#' @rdname factors
add_na_level <- function(x, name = NA, where = c("last", "first")){
  check_is_factor(x)
  where <- match.arg(where)
  lvls <- levels(x)
  if (any_na(lvls)){
    x
  } else {
    out <- unclass(x)
    n_lvls <- length(lvls)

    if (where == "first"){
      out <- out + 1L
      attr(out, "levels") <- c(name, lvls)
      out[which_na(out)] <- 1L
    } else {
      attr(out, "levels") <- c(lvls, name)
      out[which_na(out)] <- n_lvls + 1L
    }

    class(out) <- class(x)
    out
  }
}
#' @export
#' @rdname factors
drop_na_level <- function(x){

  check_is_factor(x)
  lvls <- levels(x)

  which_na_lvl <- which_na(lvls)
  if (length(which_na_lvl) == 0){
    x
  } else {
    new_lvls <- lvls[-which_na_lvl]

    matches <- collapse::fmatch(lvls, new_lvls, overid = 2L)
    out <- matches[unclass(x)]

    attributes(out) <- attributes(x)
    attr(out, "levels") <- new_lvls

    out

    ### Alternative

    # which_lvls <- seq_along(lvls)[-which_na_lvl]
    # out <- which_lvls[unclass(x)]
    #
    # attributes(out) <- attributes(x)
    # attr(out, "levels") <- levels(x)[which_lvls]
    #
    # out
  }
}

# Alternative
# drop_levels <- function(x){
#   lvls <- levels(x)
#   n_lvls <- length(lvls)
#   used_lvls <- intersect_(levels_factor(x), x)
#   if (length(used_lvls) == n_lvls){
#     x
#   } else {
#    factor_(x, levels = used_lvls)
#     # sub alternative
#     # out <- collapse::fmatch(lvls, used_lvls, overid = 2L)[unclass(x)]
#     # attributes(out) <- attributes(x)
#     # attr(out, "levels") <- as.character(used_lvls)
#   }
# }

#' @export
#' @rdname factors
drop_levels <- function(x){
  lvls <- levels(x)
  n_lvls <- length(lvls)
  which_used <- which_used_levels(x)
  if (length(which_used) == n_lvls){
    x
  } else {
    # out <- collapse::fmatch(unclass(x), which_used)
    out <- which_used[unclass(x)]
    attributes(out) <- attributes(x)
    attr(out, "levels") <- levels(x)[which_used]
    out
  }
}
#' @export
#' @rdname factors
reorder_levels <- function(x, order_by, decreasing = FALSE){
  check_is_factor(x)
  groups <- collapse::GRP(x, return.groups = FALSE)
  medians <- collapse::fmedian(order_by, g = groups, use.g.names = FALSE)

  o <- collapse::radixorderv(medians, decreasing = decreasing)
  sorted <- isTRUE(attr(o, "sorted"))
  if (sorted){
    x
  } else {
    ordered_levels <- levels(x)[o]
    factor_(x, levels = ordered_levels)
  }
  # out <- collapse::fmatch(seq_along(levels(x)), o, overid = 2L)[unclass(x)]
  # Alternative
  # out <- collapse::fmatch(x, ordered_levels, overid = 2L)

}

# rename_levels <- function(x, names){
#   check_length(names, length(levels(x)))
#   levels(x) <- names
#   x
# }


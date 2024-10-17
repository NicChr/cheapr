#' A cheaper version of `factor()` along with cheaper utilities
#'
#' @description
#' A fast version of `factor()` using the collapse package. \cr
#'
#' There are some additional utilities, most of which begin with the prefix
#' 'levels_', such as
#' `as_factor()` which is an efficient way to coerce both vectors and factors,
#' `levels_factor()` which returns the levels of a factor, as a factor,
#' `levels_used()` which returns the used levels of a factor,
#' `levels_unused()` which returns the unused levels of a factor,
#' `levels_add()` adds the specified levels onto the existing levels,
#' `levels_rm()` removes the specified levels,
#' `levels_add_na()` which adds an explicit `NA` level,
#' `levels_drop_na()` which drops the `NA` level,
#' `levels_drop()` which drops unused factor levels,
#' `levels_lump()` which returns top n levels and lumps all others into the
#' same category,
#' and finally `levels_reorder()` which reorders the levels of `x`
#' based on `y` using the ordered median values of `y` for each level.
#'
#' @returns
#' A `factor` or `character` in the case of `levels_used` and `levels_unused`.
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
#' @param n Top n number of levels to calculate.
#' @param prop Top proportion of levels to calculate.
#' This is a proportion of the total unique levels in x.
#' @param other_category Name of 'other' category.
#' @param ties Ties method to use. See `?rank`.
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
#' `levels_lump()` is a cheaper version of `forcats::lump_n()` but returns
#' levels in order of highest frequency to lowest. This can be very useful
#' for plotting.
#'
#'
#' @examples
#' library(cheapr)
#'
#' x <- factor_(sample(letters[sample.int(26, 10)], 100, TRUE), levels = letters)
#' x
#' # Used/unused levels
#'
#' levels_used(x)
#' levels_unused(x)
#'
#' # Drop unused levels
#' levels_drop(x)
#'
#' # Top 3 letters by by frequency
#' lumped_letters <- levels_lump(x, 3)
#' table(lumped_letters)
#'
#' # To remove the "other" category, use `levels_rm()`
#'
#' table(levels_rm(lumped_letters, "Other"))
#'
#' # We can use levels_lump to create a generic top n function for non-factors too
#'
#' get_top_n <- function(x, n){
#'   f <- levels_lump(factor_(x, order = FALSE), n = n)
#'   new_df(value = levels(f),
#'          count = tabulate(f))
#' }
#'
#' get_top_n(x, 3)
#'
#' # A neat way to order the levels of a factor by frequency
#' # is the following:
#'
#' table(levels_lump(x, prop = 1)) # Highest to lowest
#' table(levels_lump(x, prop = -1)) # Lowest to highest
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
  is_int64 <- inherits(x, "integer64")
  if (is_int64){
    x <- cpp_int64_to_numeric(x)
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
  } else if (is_int64){
    fct_lvls <- cpp_format_numeric_as_int64(lvls)
    } else {
    fct_lvls <- as.character(lvls)
  }
  if (inherits(x, "POSIXt") && collapse::any_duplicated(fct_lvls)){
    fct_lvls <- paste(fct_lvls, as.POSIXlt(lvls)$zone)
  }
  attr(out, "levels") <- fct_lvls
  class(out) <- c(if (ordered) "ordered" else character(), "factor")
  out
}
#' @export
#' @rdname factors
as_factor <- function(x){
  if (is.factor(x)){
    x
  } else {
    factor_(x)
  }
}

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
levels_used <- function(x){
  levels(x)[which_used_levels(x)]
}
#' @export
#' @rdname factors
levels_unused <- function(x){
  levels(x)[which_unused_levels(x)]
}
#' @export
#' @rdname factors
used_levels <- function(x){
  on.exit(.Deprecated("levels_used"))
  levels_used(x)
}
#' @export
#' @rdname factors
unused_levels <- function(x){
  on.exit(.Deprecated("levels_unused"))
  levels_unused(x)
}
#' @export
#' @rdname factors
levels_rm <- function(x, levels){
  check_is_factor(x)
  x_lvls <- levels(x)
  rm <- which_in(x_lvls, levels)

  if (length(rm) == 0){
    x
  } else {
    keep <- which_not_in(x_lvls, levels)
    factor_(x, levels = levels_factor(x)[keep])
  }
}
#' @export
#' @rdname factors
levels_add <- function(x, levels){
  check_is_factor(x)
  x_lvls <- levels(x)
  same <- which_in(levels, x_lvls)

  if (length(same) == length(levels)){
    x
  } else {
    add <- which_not_in(levels, x_lvls)
    factor_(x, levels = c(x_lvls, levels[add]))
  }
}
#' @export
#' @rdname factors
levels_add_na <- function(x, name = NA, where = c("last", "first")){
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
levels_drop_na <- function(x){
  levels_rm(x, NA)
}
#' @export
#' @rdname factors
levels_drop <- function(x){
  lvls <- levels(x)
  n_lvls <- length(lvls)
  which_used <- which_used_levels(x)
  if (length(which_used) == n_lvls){
    x
  } else {
    out <- collapse::fmatch(unclass(x), which_used, overid = 2L)
    attributes(out) <- attributes(x)
    attr(out, "levels") <- levels(x)[which_used]
    out
  }
}

#' @export
#' @rdname factors
levels_reorder <- function(x, order_by, decreasing = FALSE){
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
}
#' @export
#' @rdname factors
levels_lump <- function(x, n, prop, other_category = "Other",
                        ties = c("min", "average", "first", "last", "random", "max")){
  check_is_factor(x)
  if (!missing(n) && !missing(prop)){
    stop("Please supply either n or prop, not both")
  }
  if (!missing(prop)){
    n <- floor(prop * length(levels(x)))
  }
  check_length(n, 1)
  ties <- match.arg(ties)
  check_length(other_category, 1)
  counts <- tabulate(x, length(levels(x)))

  # Order counts
  # names(counts) <- levels(x)
  o <- order(counts, decreasing = (n >= 0))
  sorted_counts <- counts[o]
  if (ties == "min"){
    bound <- sorted_counts[min(abs(n), length(counts))]
    top <- which_(if (n >= 0) sorted_counts >= bound else sorted_counts <= bound)
  } else {
    rank <- rank(if (n >= 0) -sorted_counts else sorted_counts, ties.method = ties)
    top <- which_(rank <= abs(n))
  }
  # Does other category have > 0 items?
  other_has_counts <- length(top) != length(counts)
  if (!other_has_counts && identical(counts, sorted_counts)){
    x
  } else {
    lvls <- levels_factor(x)
    sorted_lvls <- lvls[o][top]
    # If other category already exists, we simply remove it
    # from the sorted levels
    other <- val_find(levels(sorted_lvls), other_category)
    if (length(other) > 0){
      sorted_lvls <- sorted_lvls[-val_find(sorted_lvls, other)]
      attr(sorted_lvls, "levels") <- attr(sorted_lvls, "levels")[-other]
    }
    out_levels <- factor_as_character(sorted_lvls)
    out <- collapse::fmatch(x,
                            sorted_lvls,
                            nomatch = length(sorted_lvls) + 1L,
                            overid = 2L)
    out[which_na(x)] <- NA
    if (other_has_counts){
      out_levels <- c(out_levels, other_category)
    }
    attr(out, "levels") <- out_levels
    class(out) <- "factor"
    out
  }
}
# Generic factor conversion to data representation
factor_as_type <- function(x, type){
  check_length(type, 1)

  do.call(paste0("as.", type), list(levels(x)))[unclass(x)]
}

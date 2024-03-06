#' A faster version of `factor()`
#'
#' @description
#' A fast version of `factor()` using the collapse package.
#' There are some additional utilities such as
#' `levels_factor()` which returns the levels of a factor, as a factor,
#' `used_levels()` which returns the used levels of a factor,
#' and `unused_levels()` which returns the unused levels of a factor.
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
factor_ <- function(x = integer(), levels = NULL, order = TRUE,
                    na_exclude = TRUE,
                    ordered = is.ordered(x)){
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
      lvls <- lvls[seq_len(length(lvls) - 1L)]
    } else {
      lvls <- lvls[which_not_na(lvls)]
    }
  }
  out <- collapse::fmatch(x, lvls, overid = 2L)
  fct_lvls <- as.character(lvls)
  if (inherits(x, "POSIXt") && collapse::any_duplicated(fct_lvls)){
    fct_lvls <- paste(fct_lvls, as.POSIXlt(lvls)$zone)
  }
  attr(out, "levels") <- fct_lvls
  class(out) <- c(if (ordered) "ordered" else character(),
                  # if (!na_exclude) "na_included" else character(),
                  "factor")
  out
}
#' @export
#' @rdname factors
levels_factor <- function(x){
  lvls <- levels(x)
  out <- seq_along(lvls)
  attr(out, "levels") <- lvls
  class(out) <- class(x)
  out
}
#' @export
#' @rdname factors
used_levels <- function(x){
  as.character(intersect_(levels_factor(x), x))
}
#' @export
#' @rdname factors
unused_levels <- function(x){
  as.character(setdiff_(levels_factor(x), x))
}

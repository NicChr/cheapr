#' Turn continuous data into discrete bins
#'
#' @param x A numeric vector.
#' @param breaks Break-points.
#' @param left_closed Left-closed intervals or right-closed intervals?
#' @param include_endpoint Include endpoint? Default is `FALSE`.
#' @param include_oob Include out-of-bounds values? Default is `FALSE`.
#' @param ordered Should result be an ordered factor? Default is `FALSE`.
#' @param intv_start_fun Function used to format interval start points.
#' @param intv_end_fun Function used to format interval end points.
#' @param intv_closers A length 2 character vector denoting the symbol
#' to use for closing either left or right closed intervals.
#' @param intv_openers A length 2 character vector denoting the symbol to
#' use for opening either left or right closed intervals.
#' @param intv_sep A length 1 character vector used to separate the start and
#' end points.
#' @param ... Extra arguments passed onto methods.
#'
#' @returns
#' A factor of discrete bins (intervals of start/end pairs).
#'
#' @examples
#' library(cheapr)
#'
#' # `as_discrete()` is very similar to `cut()`
#' # but more flexible as it allows you to supply
#' # formatting functions and symbols for the discrete bins
#'
#' # Here is an example of how to use the formatting functions to
#' # categorise age groups nicely
#'
#' ages <- 1:100
#'
#' age_group <- function(x, breaks){
#'   age_groups <- as_discrete(
#'     x,
#'     breaks = breaks,
#'     intv_sep = "--",
#'     intv_end_fun = function(x) x - 1,
#'     intv_openers = c("", ""),
#'     intv_closers = c("", ""),
#'     include_oob = TRUE,
#'     ordered = TRUE
#'   )
#'
#'   # Below is just renaming the last age group
#'
#'   lvls <- levels(age_groups)
#'   n_lvls <- length(lvls)
#'   max_ages <- paste0(max(breaks), "+")
#'   attr(age_groups, "levels") <- c(lvls[-n_lvls], max_ages)
#'   age_groups
#' }
#'
#' age_group(ages, seq(0, 80, 20))
#' age_group(ages, seq(0, 25, 5))
#' age_group(ages, 5)
#'
#' @rdname as_discrete
#' @export
as_discrete <- function(x, ...){
  UseMethod("as_discrete")
}
#' @rdname as_discrete
#' @export
as_discrete.numeric <- function(
    x, breaks = get_breaks(x, expand_max = TRUE),
    left_closed = TRUE,
    include_endpoint = FALSE,
    include_oob = FALSE,
    ordered = FALSE,
    intv_start_fun = prettyNum,
    intv_end_fun = prettyNum,
    intv_closers = c("[", "]"),
    intv_openers = c("(", ")"),
    intv_sep = ",",
    ...
){
  breaks <- collapse::funique(as.double(breaks), sort = TRUE)
  breaks <- na_rm(breaks)
  nb <- length(breaks)

  stopifnot(is.character(intv_closers) && length(intv_closers) == 2)
  stopifnot(is.character(intv_openers) && length(intv_openers) == 2)
  stopifnot(is.character(intv_sep) && length(intv_sep) == 1)

  # Creating labels
  if (nb < 2){
    labels <- character()
  } else {
    n <- max(nb - 1L, 0L)
    if (left_closed){
      labels <- paste0(
        intv_closers[1],
        intv_start_fun(breaks[seq_len(n)]), intv_sep,
        intv_end_fun(breaks[seq.int(to = nb, length.out = n)]),
        intv_openers[2]
      )
    } else {
      labels <- paste0(
        intv_openers[1],
        intv_end_fun(breaks[seq_len(n)]), intv_sep,
        intv_start_fun(breaks[seq.int(to = nb, length.out = n)]),
        intv_closers[2]
      )
    }
    if (anyDuplicated(labels)){
      stop("'labels' are not unique after formatting")
    }

    if (include_endpoint){
      if (left_closed && nzchar(intv_closers[2])){
        substring(labels[nb - 1L], nchar(labels[nb - 1L], "c")) <- intv_closers[2]
      } else if (nzchar(intv_closers[1])){
        substr(labels[1L], 1L, 1L) <- intv_closers[1]
      }
    }
  }

  # Binning calculation
  out <- bin(x, breaks, codes = TRUE, left_closed = left_closed,
             include_endpoint = include_endpoint,
             include_oob = include_oob)

  if (include_oob){
    if (left_closed){
      end_point <- max(breaks)
      labels <- c(labels, paste0(intv_closers[1], end_point, intv_sep, Inf, intv_openers[2]))
    } else {
      end_point <- min(breaks)
      labels <- c(paste0(intv_openers[1], -Inf, intv_sep, end_point, intv_closers[2]), labels)
    }
  }
  levels(out) <- as.character(labels)
  class(out) <- "factor"
  if (ordered) {
    class(out) <- c("ordered", class(out))
  }
  out
}
#' @export
as_discrete.integer <- as_discrete.numeric
#' @rdname as_discrete
#' @export
as_discrete.integer64 <- function(x, ...){
  as_discrete(cpp_int64_to_numeric(x), ...)
}


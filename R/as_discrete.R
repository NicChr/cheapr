#' Turn continuous data into discrete bins
#'
#' @description
#' This is a cheapr version of `cut.numeric()` which is more efficient and
#' prioritises pretty-looking breaks by default through
#' the use of `get_breaks()`.
#' Out-of-bounds values can be included naturally through the
#' `include_oob` argument. Left-closed (right-open) intervals are
#' returned by default in contrast to cut's default right-closed intervals.
#' Furthermore there is flexibility in formatting the interval bins,
#' allowing the user to specify formatting functions and symbols for
#' the interval close and open symbols.
#'
#'
#' @param x A numeric vector.
#' @param breaks Break-points.
#' @param left_closed Left-closed intervals or right-closed intervals?
#' @param include_endpoint Include endpoint? Default is `FALSE`.
#' @param include_oob Include out-of-bounds values? Default is `FALSE`.
#' This is equivalent to `breaks = c(breaks, Inf)` or
#' `breaks = c(-Inf, breaks)` when `left_closed = FALSE`.
#' If `include_endpoint = TRUE`, the endpoint interval is prioritised before
#' the out-of-bounds interval.
#' This behaviour cannot be replicated easily with `cut()`.
#' For example, these 2 expressions are not equivalent: \cr
#' \preformatted{cut(10, c(9, 10, Inf), right = F, include.lowest = T) !=
#' as_discrete(10, c(9, 10), include_endpoint = T, include_oob = T)}
#' @param ordered Should result be an ordered factor? Default is `FALSE`.
#' @param intv_start_fun Function used to format interval start points.
#' @param intv_end_fun Function used to format interval end points.
#' @param intv_closers A length 2 character vector denoting the symbol
#' to use for closing either left or right closed intervals.
#' @param intv_openers A length 2 character vector denoting the symbol to
#' use for opening either left or right closed intervals.
#' @param intv_sep A length 1 character vector used to separate the start and
#' end points.
#' @param inf_label Label to use for intervals that include infinity.
#' If left `NULL` the Unicode infinity symbol is used.
#' @param ... Extra arguments passed onto methods.
#'
#' @seealso [bin] [get_breaks]
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
#'     intv_sep = "-",
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
#' # To closely replicate `cut()` with `as_discrete()` we can use the following
#'
#' cheapr_cut <- function(x, breaks, right = TRUE,
#'                        include.lowest = FALSE,
#'                        ordered.result = FALSE){
#'   if (length(breaks) == 1){
#'     breaks <- get_breaks(x, breaks, pretty = FALSE)
#'     adj <- diff(range(breaks)) * 0.001
#'     breaks[1] <- breaks[1] - adj
#'     breaks[length(breaks)] <- breaks[length(breaks)] + adj
#'   }
#'   as_discrete(x, breaks, left_closed = !right,
#'               include_endpoint = include.lowest,
#'               ordered = ordered.result,
#'               intv_start_fun = function(x) formatC(x, digits = 3, width = 1),
#'               intv_end_fun = function(x) formatC(x, digits = 3, width = 1))
#' }
#'
#' x <- rnorm(100)
#' cheapr_cut(x, 10)
#' identical(cut(x, 10), cheapr_cut(x, 10))
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
    inf_label = NULL,
    ...
){
  breaks <- collapse::funique(as.double(breaks), sort = TRUE)
  breaks <- na_rm(breaks)
  # N breaks
  nb <- length(breaks)

  if (nb == 0){
    stop("Please provide at least 1 valid break")
  }

  # N intervals = N breaks - 1
  nintv <- max(nb - 1L, as.integer(include_endpoint))

  stopifnot(is.character(intv_closers) && length(intv_closers) == 2)
  stopifnot(is.character(intv_openers) && length(intv_openers) == 2)
  stopifnot(is.character(intv_sep) && length(intv_sep) == 1)

  # Creating labels
  if (nb < (2 - include_endpoint)){
    labels <- character()
  } else {

    if (left_closed){
      labels <- paste0(
        intv_closers[1],
        intv_start_fun(breaks[seq_len(nintv)]), intv_sep,
        intv_end_fun(breaks[seq.int(to = nb, length.out = nintv)]),
        intv_openers[2]
      )
    } else {
      labels <- paste0(
        intv_openers[1],
        intv_end_fun(breaks[seq_len(nintv)]), intv_sep,
        intv_start_fun(breaks[seq.int(to = nb, length.out = nintv)]),
        intv_closers[2]
      )
    }
    if (anyDuplicated(labels)){
      stop("'labels' are not unique after formatting")
    }

    if (include_endpoint && nb >= 1){
      if (left_closed && nzchar(intv_closers[2])){
        substring(labels[nintv], nchar(labels[nintv], "c")) <- intv_closers[2]
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
    if (is.null(inf_label)){
      inf_label <- "\u221E"
    }
    if (left_closed){
      end_point <- max(breaks)
      labels <- c(labels, paste0(intv_closers[1], end_point, intv_sep, inf_label, intv_openers[2]))
    } else {
      end_point <- min(breaks)
      labels <- c(paste0(intv_openers[1], "-", inf_label, intv_sep, end_point, intv_closers[2]), labels)
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


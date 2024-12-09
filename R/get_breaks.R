#' Pretty break-points for continuous (numeric) data
#'
#' @description
#' The distances between break-points are always equal in this implementation.
#'
#' @param x A numeric vector.
#' @param n Number of breakpoints. You may get less or more than requested.
#' @param pretty Should pretty break-points be prioritised? Default is `TRUE`.
#' If `FALSE` bin-widths will be calculated as `diff(range(x)) / n`.
#' @param expand_min Should smallest break be extended beyond the
#' minimum of the data? Default is `FALSE`.
#' If `TRUE` then `min(get_breaks(x))` is ensured to be less than `min(x)`.
#' @param expand_max Should largest break be extended beyond the maximum
#' of the data? Default is `TRUE`.
#' If `TRUE` then `max(get_breaks(x))` is ensured
#' to be greater than `max(x)`.
#' @param ... Extra arguments passed onto methods.
#'
#' @seealso [bin] [as_discrete]
#'
#' @returns
#' A numeric vector of break-points.
#'
#' @examples
#' library(cheapr)
#'
#' set.seed(123)
#' ages <- sample(0:80, 100, TRUE)
#'
#' # Pretty
#' get_breaks(ages, n = 10)
#' # Not-pretty
#' # bin-width is diff(range(ages)) / n_breaks
#' get_breaks(ages, n = 10, pretty = FALSE)
#'
#' # `get_breaks()` is left-biased in a sense, meaning that
#' # the first break is always <= `min(x)` but the last break
#' # may be < `max(x)`
#'
#' # To get right-biased breaks we can use a helper like so..
#'
#' right_breaks <- function(x, ...){
#'   -get_breaks(-x, ...)
#' }
#'
#' get_breaks(4:24, 10)
#' right_breaks(4:24, 10)
#'
#' # Use `rev()` to ensure they are in ascending order
#' rev(right_breaks(4:24, 10))
#'
#' @rdname get_breaks
#' @export
get_breaks <- function(x, n = 10, ...){
  UseMethod("get_breaks")
}
#' @rdname get_breaks
#' @export
get_breaks.numeric <- function(x, n = 10,
                               pretty = TRUE,
                               expand_min = FALSE,
                               expand_max = pretty,
                               ...){
  rng <- as.double(collapse::frange(x, na.rm = TRUE, finite = TRUE))
  cpp_fixed_width_breaks(
    rng[1], rng[2], as.double(n),
    as.logical(pretty), as.logical(expand_min), as.logical(expand_max)
  )
}

#' @export
get_breaks.integer <- get_breaks.numeric
#' @rdname get_breaks
#' @export
get_breaks.integer64 <- function(x, n = 10, ...){
  get_breaks(cpp_int64_to_numeric(x), n, ...)
}

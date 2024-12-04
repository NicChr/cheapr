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
                               expand_max = TRUE,
                               ...){
  check_length(n, 1L)
  stopifnot(n >= 1)
  stopifnot(is.finite(n))

  rng <- as.double(collapse::frange(x, na.rm = TRUE, finite = TRUE))
  if (any_na(rng) || all(is.infinite(rng))){
    return(NA_real_)
  }
  start <- rng[1]
  end <- rng[2]
  rng_width <- end - start
  spans_zero <- abs(diff(sign(rng))) == 2
  zero_range <- isTRUE(rng_width == 0)
  tol <- sqrt(.Machine$double.eps)

  if (zero_range){
    if (isTRUE(start == 0)){
      rng_width <- 1
    } else {
      rng_width <- abs(start)
    }
    return(
      seq(start - (rng_width / 1000),
          end + (rng_width / 1000),
          length.out = (n + 1))
    )
  }

  if (!pretty){
    width <- rng_width / n
    if (expand_min){
      start <- start - width
    }
    if (expand_max){
      end <- end + width
    }
    seq(start, end, by = width)

  } else {

    # Making breaks prettier

    # This is the number of orders of magnitude the data spans
    scale_diff <- log_scale(rng_width)

    # If large range & relatively small starting value
    # floor start to the nearest difference in orders of magnitude

    # If range spans across zero and start val is small
    if (scale_diff >= 1 && spans_zero && abs(start) < 1){
      adj_start <- nearest_floor(start, 10^(log_scale(end)))
    } else {
      adj_start <- nearest_floor(start, 10^(scale_diff))
    }

    # Calculate bin-width (guaranteed to span end-points inclusively)
    adj_rng_width <- end - adj_start
    bin_width <- adj_rng_width / n

    # Make width look nicer
    if (bin_width > 2 && bin_width < 5){
      adj_width <- 5
    } else if (bin_width > 5 && bin_width < 10){
      adj_width <- 10
    } else if (bin_width < 1){
      adj_width <- nearest_ceiling(bin_width, (10^(ceiling(log10(bin_width)))) / 2)
      } else {
      adj_width <- pretty_ceiling(bin_width)
      }

    # n is our first best guess at number of final breaks
    n_breaks <- n


    # These functions only work on length-1 vectors

    # Is x a whole number?
    is_whole_number <- function(x, .tol = tol){
      if (length(x) != 1){
        stop("x must be a length-1 vector")
      }
      isTRUE(abs(x - round(x)) < .tol)
    }

    # Integer division but if division is
    # very close to a whole number then we round
    approx_int_div <- function(x, y, .tol = tol){
      out <- x / y
      if (is_whole_number(out, .tol)) round(out) else trunc(out)
    }

    # If breaks span zero make sure they actually land on zero
    if (spans_zero){
      zero_adjustment <- adj_start + (adj_width * ceiling(abs(adj_start) / adj_width))
      lands_on_zero <- abs(zero_adjustment) < tol
      if (!lands_on_zero){
        adj_start <- adj_start - zero_adjustment
        n_breaks <- n_breaks + approx_int_div(zero_adjustment, adj_width)
      }
    }

    # Final break?
    adj_end <- seq_to(n_breaks, adj_start, by = adj_width)

    # If too many breaks, reduce
    if (adj_end > end){
      n_rm <- ceiling( (adj_end - end) / adj_width)
      n_breaks <- n_breaks - n_rm
      adj_end <- adj_end - (adj_width * n_rm)
    }

    # adj_end might also have too few breaks
    if (adj_end < end){
      n_add <- approx_int_div(end - adj_end, adj_width)
      n_breaks <- n_breaks + n_add
      adj_end <- adj_end + (adj_width * n_add)
    }
    if (adj_start < start){
      n_rm <- approx_int_div(start - adj_start, adj_width)
      n_breaks <- n_breaks - n_rm
      adj_start <- adj_start + (adj_width * n_rm)
    }

    # At this point, adj_end will be in range (end - adj_width, end]
    # If we want last break to extend beyond data, add adj_width to it
    if (expand_max && adj_end <= end){
      n_breaks <- n_breaks + 1
      adj_end <- adj_end + adj_width
    }

    # Reduce floating-point error
    # If width is almost a whole number, just round it

    if (is_whole_number(adj_width)){
      adj_width <- as.integer(round(adj_width))
    }
    if (is_whole_number(adj_start)){
      adj_start <- as.integer(round(adj_start))
    }

    # If last break is >= 2^31 use doubles
    if (!isTRUE(is_integerable(adj_end))){
      adj_start <- as.double(adj_start)
      adj_width <- as.double(adj_width)
    }

    if (expand_min && adj_start >= start){
      n_breaks <- n_breaks + 1
      adj_start <- adj_start - adj_width
    }
    seq(adj_start, by = adj_width, length.out = n_breaks)
  }
}

#' @export
get_breaks.integer <- get_breaks.numeric
#' @rdname get_breaks
#' @export
get_breaks.integer64 <- function(x, n = 10, ...){
  get_breaks(cpp_int64_to_numeric(x), n, ...)
}

log_scale <- function(x, base = 10){
  y <- val_replace(x, 0, 1)
  floor(log(abs(y), base = base))
}
nearest_floor <- function(x, n){
  floor(x / n) * n
}
nearest_ceiling <- function(x, n){
  ceiling(x / n) * n
}

## Based on the log10 scale of x,
## it is rounded down/up to the nearest order of magnitude

pretty_floor <- function(x, base = 10){
  scale <- log_scale(x, base)
  nearest_floor(x, (base ^ scale))
}
pretty_ceiling <- function(x, base = 10){
  scale <- log_scale(x, base)
  nearest_ceiling(x, (base ^ scale))
}

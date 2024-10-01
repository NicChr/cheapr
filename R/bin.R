#' A sometimes cheaper but argument richer alternative to `.bincode()`
#'
#' @description
#' When `x` is an integer vector, `bin()` is cheaper than `.bincode()` as
#' no coercion to a double vector occurs. This alternative also has more
#' arguments that allow you to return the start values of the binned vector,
#' as well as including out-of-bounds intervals.
#'
#' @param x A numeric vector.
#' @param breaks A numeric vector of breaks.
#' @param right Should intervals be right-closed (and left-open)?
#' Default is `FALSE` which means they are right-open (and left-closed).
#' @param include_lowest See `?.bincode` for details.
#' @param include_oob Should out-of-bounds interval be included?
#' Default is `FALSE`. This is the equivalent of adding `Inf` as
#' the last value of the breaks, or the
#' first value of the breaks if `right = TRUE`.
#' @param codes Should an integer vector indicating which bin the values
#' fall into be returned? Default is `TRUE`. If `FALSE` the
#' start values of the respective bin intervals are returned, i.e the
#' corresponding breaks.
#'
#' @returns
#' Either an integer vector of codes indicating which bin the values fall
#' into, or the start of the intervals for which each value falls into.
#'
#' @export
bin <- function(x, breaks, right = FALSE, include_lowest = FALSE,
                include_oob = FALSE, codes = TRUE){
  cpp_bin(x, breaks, right = right, include_lowest = include_lowest,
          include_oob = include_oob, codes = codes)
}

## Based on the log10 scale of x, it is rounded down/up to the nearest order
## of magnitude

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
pretty_floor <- function(x, base = 10){
  scale <- log_scale(x, base)
  nearest_floor(x, (base ^ scale))
}
pretty_ceiling <- function(x, base = 10){
  scale <- log_scale(x, base)
  nearest_ceiling(x, (base ^ scale))
}

# get_breaks <- function(x, n = 5, pretty = TRUE, extend_range = FALSE){
#   check_length(n, 1L)
#   stopifnot(n >= 1)
#
#   if (!is.numeric(x)){
#     stop("'x' must be numeric")
#   }
#
#   rng <- as.double(collapse::frange(x, na.rm = TRUE))
#   rng_width <- diff(rng)
#   start <- rng[1]
#   end <- rng[2]
#
#   if (rng_width == 0){
#     return(seq(start - 0.05, end + 0.05, by = 0.1 / max((n - 1), 1)))
#   }
#
#   if (!pretty){
#
#     # The way cut() does it
#     width <- rng_width / n
#     out <- seq(start, end, by = width)
#     if (extend_range){
#       adj_start <- start - (rng_width * 0.001)
#       adj_end <- end + (rng_width * 0.001)
#       out[1] <- adj_start
#       out[length(out)] <- adj_end
#     }
#     out
#
#   } else {
#
#     # Make end-points prettier
#
#     adj_start <- round(pretty_floor(start), 6)
#     adj_rng_width <- end - adj_start
#
#     # Calculate bin-width (guaranteed to span end-points inclusively)
#     bin_width <- adj_rng_width / n
#     adj_width <- round(pretty_ceiling(bin_width), 6)
#
#     # Reduce floating-point error
#     # If width is almost a whole number, just round it
#
#     if (all(abs(adj_width - round(adj_width)) < sqrt(.Machine$double.eps), na.rm = TRUE)) {
#       adj_width <- round(adj_width)
#     }
#
#     # Breaks
#
#     n_breaks <- n + 1
#     if (extend_range){
#       n_breaks <- n_breaks + 2 # 2 for good measure
#     }
#
#     if (extend_range && adj_start >= start){
#       adj_start <- adj_start - adj_width
#     }
#
#     out <- seq(adj_start, by = adj_width, length.out = n_breaks)
#
#     rm <- which_(out > end)
#
#     if (extend_range && (length(rm) > 0)){
#       rm <- rm[-1L]
#     }
#
#     # Remove outlier breaks
#
#     if (length(rm) > 0){
#       out <- out[-rm]
#     }
#
#     out
#   }
# }

get_breaks <- function(x, n = 5, pretty = TRUE,
                       expand_min = FALSE,
                       expand_max = FALSE){
  check_length(n, 1L)
  stopifnot(n >= 1)

  if (!is.numeric(x)){
    stop("'x' must be numeric")
  }

  rng <- as.double(collapse::frange(x, na.rm = TRUE))
  rng_width <- diff(rng)
  start <- rng[1]
  end <- rng[2]

  if (rng_width == 0){
    return(seq(start - 0.05, end + 0.05, by = 0.1 / max((n - 1), 1)))
  }

  if (!pretty){

    # The way cut() does it
    width <- rng_width / n
    out <- seq(start, end, by = width)
    if (expand_min){
      adj_start <- start - (rng_width * 0.001)
      out[1] <- adj_start
    }
    if (expand_max){
      adj_end <- end + (rng_width * 0.001)
      out[length(out)] <- adj_end
    }
    out

  } else {

    # Make end-points prettier

    adj_start <- round(pretty_floor(start), 6)
    adj_rng_width <- end - adj_start

    # Calculate bin-width (guaranteed to span end-points inclusively)
    bin_width <- adj_rng_width / n
    adj_width <- round(pretty_ceiling(bin_width), 6)

    # Reduce floating-point error
    # If width is almost a whole number, just round it

    if (all(abs(adj_width - round(adj_width)) < sqrt(.Machine$double.eps), na.rm = TRUE)) {
      adj_width <- round(adj_width)
    }

    # Breaks

    n_breaks <- n + 1
    if (expand_max){
      n_breaks <- n_breaks + 2 # 2 for good measure
    }

    if (expand_min && adj_start >= start){
      adj_start <- adj_start - adj_width
    }

    out <- seq(adj_start, by = adj_width, length.out = n_breaks)

    rm <- which_(out > end)

    if (expand_max && (length(rm) > 0)){
      rm <- rm[-1L]
    }

    # Remove outlier breaks

    if (length(rm) > 0){
      out <- out[-rm]
    }

    out
  }
}
make_discrete <- function(
    x, breaks = get_breaks(x, 7, expand_max = TRUE),
    labels_fun = prettyNum,
    left_closed = TRUE,
    include_lowest = FALSE,
    include_oob = FALSE
){
  if (!is.numeric(x)){
    stop("'x' must be numeric")
  }
  right <- !left_closed
  breaks <- sort.int(as.double(breaks))
  nb <- length(breaks)
  # Creating labels
  if (nb < 2){
    labels <- character()
  } else {
    n <- max(nb - 1L, 0L)
    formatted_breaks <- match.fun(labels_fun)(breaks)
    if (length(formatted_breaks) != length(breaks)){
      stop("`labels_fun` must produce a vector of length equal to the length of `breaks`")
    }
    if (anyDuplicated(formatted_breaks)){
      stop("'breaks' are not unique either before or after formatting")
    }
    if (left_closed){
      labels <- paste0(
        "[",
        formatted_breaks[seq_len(n)], ",",
        formatted_breaks[seq.int(to = nb, length.out = n)],
        ")"
      )
    } else {
      labels <- paste0(
        "(",
        formatted_breaks[seq_len(n)], ",",
        formatted_breaks[seq.int(to = nb, length.out = n)],
        "]"
      )
    }

    if (include_lowest && all(formatted_breaks[-1L] != formatted_breaks[-nb])){
      if (left_closed){
        substring(labels[nb - 1L], nchar(labels[nb - 1L], "c")) <- "]"
      } else {
        substr(labels[1L], 1L, 1L) <- "["
      }
    }
  }
  out <- bin(x, breaks, codes = TRUE, right = right,
             include_lowest = include_lowest,
             include_oob = include_oob)
  if (include_oob){
    if (left_closed){
      end_point <- max(breaks)
      if ( (include_lowest && any(x > end_point)) ||
           (!include_lowest && any(x >= end_point) )){
        labels <- c(labels, paste0("[", end_point, ",", Inf, ")"))
      }
    } else {
      end_point <- min(breaks)
      if ( (include_lowest && any(x < end_point, TRUE) ) ||
           (!include_lowest && any(x <= end_point) )){
        labels <- c(paste0("(", -Inf, ",", end_point, "]"), labels)
      }
    }
  }
  levels(out) <- as.character(labels)
  class(out) <- "factor"
  out
}

# get_breaks <- function(x, n = 5, pretty = TRUE){
#   check_length(n, 1L)
#   stopifnot(n >= 1)
#
#   if (!is.numeric(x)){
#     stop("'x' must be numeric")
#   }
#
#   rng <- as.double(collapse::frange(x, na.rm = TRUE))
#   rng_width <- diff(rng)
#   start <- rng[1]
#   end <- rng[2]
#
#   if (rng_width == 0){
#     return(seq(start - 0.05, end + 0.05, by = 0.1 / max((n - 1), 1)))
#   }
#
#   if (!pretty){
#
#     # The way cut() does it
#     width <- rng_width / n
#     out <- seq(start, end, by = width)
#     adj_start <- start - (rng_width * 0.001)
#     adj_end <- end + (rng_width * 0.001)
#     out[1] <- adj_start
#     out[length(out)] <- adj_end
#     out
#
#   } else {
#
#     # Make end-points prettier
#
#     adj_start <- round(pretty_floor(start), 6)
#     adj_end <- round(pretty_ceiling(end), 6)
#     adj_rng_width <- adj_end - adj_start
#
#     # Calculate bin-width (guaranteed to span end-points inclusively)
#     bin_width <- adj_rng_width / n
#     adj_width <- round(pretty_ceiling(bin_width), 6)
#
#     # Reduce floating-point error
#     # If width is almost a whole number, just round it
#
#     if (all(abs(adj_width - round(adj_width)) < sqrt(.Machine$double.eps), na.rm = TRUE)) {
#       adj_width <- round(adj_width)
#     }
#
#     # We want breaks to span entire data, most likely going out-of-bounds
#     n_breaks <- seq_size(adj_start, adj_end, by = adj_width)
#     if (length(n_breaks) == 0) {
#       n_breaks <- 0
#     }
#     last_break <- seq_to(n_breaks, adj_start, adj_width)
#     if (last_break < adj_end){
#       n_breaks <- n_breaks + 1
#     }
#
#     # Breaks
#     out <- seq(adj_start, by = adj_width, length.out = n_breaks)
#
#     # Remove outlier breaks
#
#     # We are fine with at most 1 break > max(x)
#
#     outlier_loc <- which_(out > end)
#     rm <- outlier_loc[-1L]
#     if (length(rm) > 0){
#       out <- out[-rm]
#     }
#
#     out
#   }
# }
# get_breaks2 <- function(x, n = 5, pretty = TRUE, extend_range = FALSE){
#   check_length(n, 1L)
#   stopifnot(n >= 1)
#
#   if (!is.numeric(x)){
#     stop("'x' must be numeric")
#   }
#
#   rng <- as.double(collapse::frange(x, na.rm = TRUE))
#   rng_width <- diff(rng)
#   start <- rng[1]
#   end <- rng[2]
#
#   if (rng_width == 0){
#     return(seq(start - 0.05, end + 0.05, by = 0.1 / max((n - 1), 1)))
#   }
#
#   if (!pretty){
#
#     # The way cut() does it
#     width <- rng_width / n
#     out <- seq(start, end, by = width)
#     if (extend_range){
#       adj_start <- start - (rng_width * 0.001)
#       adj_end <- end + (rng_width * 0.001)
#       out[1] <- adj_start
#       out[length(out)] <- adj_end
#     }
#     out
#
#   } else {
#
#     # Make end-points prettier
#
#     adj_start <- round(pretty_floor(start), 6)
#     adj_rng_width <- end - adj_start
#
#     # Calculate bin-width (guaranteed to span end-points inclusively)
#     bin_width <- adj_rng_width / n
#     adj_width <- round(pretty_ceiling(bin_width), 6)
#
#     # If abs(adj_width) is in interval [1, 10)
#     # We might want to increment using a number to the nearest 5
#
#     # if (log_scale(adj_width) == 0){
#     #   adj_width <- nearest_ceiling(adj_width, 5)
#     # }
#
#     # Reduce floating-point error
#     # If width is almost a whole number, just round it
#
#     if (all(abs(adj_width - round(adj_width)) < sqrt(.Machine$double.eps), na.rm = TRUE)) {
#       adj_width <- round(adj_width)
#     }
#
#     # We want breaks to span entire data, most likely going out-of-bounds
#     # adj_end <- seq_to(n, adj_start, adj_width)
#     #
#     # n_breaks <- n
#     # if (adj_end < end){
#     #   n_breaks <- n_breaks + 1
#     # }
#
#     # Breaks
#
#     # Generally better to start breaks at 0 if min(x) == 1
#     # if (adj_start == 1){
#     #   adj_start <- 0
#     # }
#     n_breaks <- n + 1
#     if (extend_range){
#      n_breaks <- n_breaks + 2 # 2 for good measure
#     }
#
#     if (extend_range && adj_start >= start){
#       adj_start <- adj_start - adj_width
#     }
#
#     out <- seq(adj_start, by = adj_width, length.out = n_breaks)
#
#     rm <- which_(out > end)
#
#     if (extend_range && (length(rm) > 0)){
#       rm <- rm[-1L]
#     }
#
#     # Remove outlier breaks
#
#     if (length(rm) > 0){
#       out <- out[-rm]
#     }
#
#     out
#   }
# }

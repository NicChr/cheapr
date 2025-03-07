#' An alternative to `summary()` inspired by the skimr package
#'
#' @description
#' A cheaper `summary()` function, designed for larger data.
#'
#' @param x A vector or data frame.
#' @param hist Should in-line histograms be returned? Default is `FALSE`.
#' @param digits How many decimal places should the summary statistics be
#' printed as? Default is 2.
#' @param ... Further arguments passed onto methods. Currently unused.
#'
#' @returns
#' An object of class "overview".
#' Under the hood this is just a list of data frames.
#' Key summary statistics are reported in each data frame.
#'
#' @details
#' No rounding of statistics is done except in printing which can be controlled
#' either through the `digits` argument in `overview()`, or by setting the
#' option `options(cheapr.digits)`. \cr
#' To access the underlying data, for example the numeric summary,
#' just use `$numeric`, e.g. `overview(rnorm(30))$numeric`.
#'
#' @examples
#' library(cheapr)
#' overview(iris)
#'
#' # With histograms
#' overview(airquality, hist = TRUE)
#'
#' # Round to 0 decimal places
#' overview(airquality, digits = 0)
#'
#' # We can set an option for all overviews
#' options(cheapr.digits = 1)
#' overview(rnorm(100))
#' options(cheapr.digits = 2) # The default
#'
#' @rdname overview
#' @export
overview <- function(x, digits = getOption("cheapr.digits", 2), ...){
  UseMethod("overview")
}
#' @rdname overview
#' @export
overview.default <- function(x, digits = getOption("cheapr.digits", 2), ...){
  overview(new_df(x = x), digits = digits, ...)
}
#' @rdname overview
#' @export
overview.logical <- function(x, digits = getOption("cheapr.digits", 2), ...){
  overview(new_df(x = as.logical(x)), digits = digits, ...)
}
#' @rdname overview
#' @export
overview.integer <- function(x, digits = getOption("cheapr.digits", 2), hist = TRUE, ...){
  out <- overview(new_df(x = as.integer(x)), hist = hist, digits = digits, ...)
  out$cols <- NA_integer_
  out
}
#' @rdname overview
#' @export
overview.numeric <- function(x, digits = getOption("cheapr.digits", 2), hist = TRUE, ...){
  out <- overview(new_df(x = as.numeric(x)), hist = hist, digits = digits, ...)
  out$cols <- NA_integer_
  out
}
#' @rdname overview
#' @export
overview.integer64 <- function(x, digits = getOption("cheapr.digits", 2), hist = TRUE, ...){
  out <- overview(cpp_int64_to_numeric(x), hist = hist, digits = digits, ...)
  out$numeric$class <- class(x)[1]
  out
}
#' @rdname overview
#' @export
overview.character <- function(x, digits = getOption("cheapr.digits", 2), ...){
  out <- overview(new_df(x = as.character(x)), digits = digits, ...)
  out$cols <- NA_integer_
  out
}
#' @rdname overview
#' @export
overview.factor <- function(x, digits = getOption("cheapr.digits", 2), ...){
  out <- overview(new_df(x = as.factor(x)), digits = digits, ...)
  out$cols <- NA_integer_
  out
}
#' @rdname overview
#' @export
overview.Date <- function(x, digits = getOption("cheapr.digits", 2), ...){
  out <- overview(new_df(x = as.Date(x)), digits = digits, ...)
  out$cols <- NA_integer_
  out
}
#' @rdname overview
#' @export
overview.POSIXt <- function(x, digits = getOption("cheapr.digits", 2), ...){
  out <- overview(new_df(x = as.POSIXct(x)), digits = digits, ...)
  out$cols <- NA_integer_
  out
}
#' @rdname overview
#' @export
overview.ts <- function(x, digits = getOption("cheapr.digits", 2), ...){
  ts_df <- as.data.frame(x)
  out <- overview(transform_all(ts_df, as.numeric), digits = digits, ...)
  out$time_series <- out$numeric
  out$numeric <- sset(out$numeric, 0)
  out$time_series$class <- class(x)[1]
  out$time_series$growth_rate <- pluck_row(summarise_all(ts_df, cpp_growth_rate))
  out
}
#' @rdname overview
#' @export
overview.zoo <- overview.ts
#' @rdname overview
#' @export
overview.data.frame <- function(x, digits = getOption("cheapr.digits", 2), hist = TRUE, ...){
  N <- nrow(x)
  num_cols <- ncol(x)
  skim_df <- x
  data_nms <- names(skim_df)

  col_classes <- vapply(skim_df, function(x) class(x)[length(class(x))], "")
  out <- as_df(cheapr_enframe(col_classes, name = "col", value = "class"))
  chr_vars <- data_nms[vapply(skim_df, is.character, FALSE,
                              USE.NAMES = FALSE)]
  skim_df <- transform_all(
    skim_df, function(x) factor_(x, ordered = TRUE), .cols = chr_vars
  )
  lgl_vars <- data_nms[vapply(skim_df, is.logical, FALSE)]
  num_vars <- data_nms[vapply(skim_df, function(x) inherits(x,
                                                            c("integer", "numeric", "integer64")), FALSE)]
  ### SUBSET OF NUM_VARS
  int64_vars <- data_nms[vapply(skim_df, function(x) inherits(x, "integer64"), FALSE)]

  date_vars <- data_nms[vapply(skim_df, function(x) inherits(x, "Date"), FALSE)]
  datetime_vars <- data_nms[vapply(skim_df, function(x) inherits(x, "POSIXt"),  FALSE)]
  ts_vars <- data_nms[vapply(skim_df, function(x) inherits(x, c("ts", "zoo")),  FALSE)]
  cat_vars <- data_nms[vapply(skim_df, is.factor, FALSE)]
  other_vars <- setdiff(data_nms, c(lgl_vars, num_vars, date_vars,
                                    datetime_vars, ts_vars, cat_vars))

# Logical -----------------------------------------------------------------

  lgl_data <- sset_col(skim_df, lgl_vars)
  which_lgl <- which_in(out[["col"]], lgl_vars)
  lgl_out <- sset(out, which_lgl)
  lgl_out <- df_add_cols(
    lgl_out,
    list(
      n_missing = NA_integer_,
      p_complete = NA_real_,
      n_true = NA_integer_,
      n_false = NA_integer_,
      p_true = NA_real_
    )
  )
  if (N > 0L && length(which_lgl) > 0) {
    lgl_out$n_missing <- pluck_row(summarise_all(lgl_data, na_count), 1)
    lgl_out$p_complete <- pluck_row(summarise_all(lgl_data, prop_complete), 1)
    lgl_out$n_true <- pluck_row(summarise_all(lgl_data, function(x) sum(x, na.rm = TRUE)), 1)
    lgl_out$n_false <- N - lgl_out[["n_missing"]] - lgl_out[["n_true"]]
    lgl_out$p_true <- lgl_out[["n_true"]] / (N - lgl_out[["n_missing"]])
  }

# Numeric -----------------------------------------------------------------

  num_data <- sset_col(skim_df, num_vars)

  ## Coerce int64 to double
  num_data <- transform_all(num_data, cpp_int64_to_numeric, int64_vars)

  which_num <- which_in(out[["col"]], num_vars)
  num_out <- sset(out, which_num)
  num_out <- df_add_cols(
    num_out,
    list(
      n_missing = NA_integer_,
      p_complete = NA_real_,
      n_unique = NA_integer_,
      mean = NA_real_,
      p0 = NA_real_,
      p25 = NA_real_,
      p50 = NA_real_,
      p75 = NA_real_,
      p100 = NA_real_,
      iqr = NA_real_,
      sd = NA_real_
    )
  )

  if (hist){
    num_out$hist <- character(length(which_num))
  }

  if (N > 0L && length(which_num) > 0) {
    num_summaries <- unname(lapply(num_data, numeric_summary))
    num_out <- df_add_cols(
      num_out, collapse::rowbind(num_summaries)
    )
    if (hist){
      num_out$hist <- pluck_row(summarise_all(
        num_data, inline_hist
      ), 1)
    }
  }

# Dates -------------------------------------------------------------------

  date_data <- sset_col(skim_df, date_vars)
  which_date <- which_in(out[["col"]], date_vars)
  date_out <- sset(out, which_date)
  date_out <- df_add_cols(
    date_out,
    list(
      n_missing = NA_integer_,
      p_complete = NA_real_,
      n_unique = NA_integer_,
      min = NA_real_,
      max = NA_real_
    )
  )
  if (N > 0L && length(which_date) > 0) {
    date_out$n_missing <- pluck_row(summarise_all(date_data, na_count), 1)
    date_out$p_complete <- pluck_row(summarise_all(date_data, prop_complete), 1)
    date_out$n_unique <- pluck_row(summarise_all(date_data, n_unique), 1)
    date_out$n_unique <- date_out$n_unique - (date_out$n_missing > 0L)
    date_out$min <- pluck_row(summarise_all(date_data, collapse::fmin), 1)
    date_out$max <- pluck_row(summarise_all(date_data, collapse::fmax), 1)
    class(date_out$min) <- "Date"
    class(date_out$max) <- "Date"
  }


# Date-Times --------------------------------------------------------------

  datetime_data <- sset_col(skim_df, datetime_vars)
  datetime_data <- transform_all(datetime_data, as.POSIXct)
  which_datetime <- which_in(out[["col"]], datetime_vars)
  datetime_out <- sset(out, which_datetime)
  datetime_out <- df_add_cols(
    datetime_out,
    list(
      n_missing = NA_integer_,
      p_complete = NA_real_,
      n_unique = NA_integer_,
      tzone = NA_character_,
      min = NA_real_,
      max = NA_real_
    )
  )
  if (N > 0L && length(which_datetime) > 0) {
    datetime_out$n_missing <- pluck_row(summarise_all(datetime_data, na_count), 1)
    datetime_out$p_complete <- pluck_row(summarise_all(datetime_data, prop_complete), 1)
    datetime_out$n_unique <- pluck_row(summarise_all(datetime_data, n_unique), 1)
    datetime_out$n_unique <- datetime_out$n_unique - (datetime_out$n_missing > 0L)
    datetime_out$min <- pluck_row(summarise_all(datetime_data, collapse::fmin), 1)
    datetime_out$max <- pluck_row(summarise_all(datetime_data, collapse::fmax), 1)
    datetime_out$min <- .POSIXct(datetime_out$min, tz = "UTC")
    datetime_out$max <- .POSIXct(datetime_out$max, tz = "UTC")
    for (i in seq_len(nrow(datetime_out))){
      datetime_out[["tzone"]][i] <- tzone(skim_df[[datetime_out[["col"]][i]]])
    }
  }

  # Time-Series -----------------------------------------------------------------

  ts_data <- sset_col(skim_df, ts_vars)
  which_ts <- which_in(out[["col"]], ts_vars)
  ts_out <- sset(out, which_ts)
  if (N > 0L && length(which_ts) > 0) {
    ts_overviews <- new_list(nrow(ts_out))
    for (i in seq_along(ts_overviews)){
      ts_overviews[[i]] <- overview(
        ts_data[[ts_out[["col"]][i]]],
        hist = hist, digits = digits, ...
      )$time_series
      if (length(attr(ts_overviews[[i]], "row.names")) > 1){
        ts_overviews[[i]][["col"]] <- paste0(ts_out[["col"]][i], "_",
                                             ts_overviews[[i]][["col"]])
      } else {
        ts_overviews[[i]][["col"]] <- ts_vars[i]
      }
    }
    ts_out <- collapse::rowbind(ts_overviews)
  }

# Categorical -------------------------------------------------------------

  cat_data <- sset_col(skim_df, cat_vars)
  which_cat <- which_in(out[["col"]], cat_vars)
  cat_out <- sset(out, which_cat)
  cat_out <- df_add_cols(
    cat_out,
    list(
      n_missing = NA_integer_,
      p_complete = NA_real_,
      n_unique = NA_integer_,
      n_levels = NA_integer_,
      min = NA_character_,
      max = NA_character_
    )
  )
  if (N > 0L && length(which_cat) > 0) {
    cat_out$n_missing <- pluck_row(summarise_all(cat_data, na_count), 1)
    cat_out$p_complete <- pluck_row(summarise_all(cat_data, prop_complete), 1)
    cat_out$n_unique <- pluck_row(summarise_all(cat_data, n_unique), 1)
    cat_out$n_unique <- cat_out$n_unique - (cat_out$n_missing > 0L)
    cat_out$min <- pluck_row(summarise_all(cat_data, collapse::fmin), 1)
    cat_out$max <- pluck_row(summarise_all(cat_data, collapse::fmax), 1)
    cat_out$min <- as.character(cat_out$min)
    cat_out$max <- as.character(cat_out$max)
    for (i in seq_len(nrow(cat_out))){
      if (cat_out[["class"]][i] == "factor"){
        cat_out[["n_levels"]][i] <- length(levels(skim_df[[cat_out[["col"]][i]]]))
      }
    }
  }

  # Other -------------------------------------------------------------

  other_data <- sset_col(skim_df, other_vars)
  which_other <- which_in(out[["col"]], other_vars)
  other_out <- sset(out, which_other)
  other_out <- df_add_cols(
    other_out,
    list(
      n_missing = NA_integer_,
      p_complete = NA_real_,
      n_unique = NA_integer_
    )
  )
  if (N > 0L && length(which_other) > 0) {
    other_out$n_missing <- pluck_row(summarise_all(
      other_data, function(x) na_count(x, recursive = FALSE)
      ), 1)
    other_out$p_complete <- pluck_row(summarise_all(
      other_data, function(x) prop_complete(x, recursive = FALSE)
      ), 1)
    other_out$n_unique <- pluck_row(summarise_all(
      other_data, function(x) length(unique(x))
      ), 1)
    other_out$n_unique <- other_out$n_unique - (other_out$n_missing > 0L)
  }

  out <- list(
    obs = N, cols = num_cols,
    print_digits = digits,
    logical = lgl_out,
    numeric = num_out,
    date = date_out,
    datetime = datetime_out,
    time_series = ts_out,
    categorical = cat_out,
    other = other_out
  )
  class(out) <- "overview"
  out
}
#' @export
print.overview <- function(x, max = NULL, ...){
  temp <- unclass(x)
  digits <- temp[["print_digits"]]

  # A nice rule for rounding prettily
  # For numbers in { |x| < 1 }
  # use significant digits
  # Otherwise use decimal digits

  pretty_num <- function(
    x, round_digits = digits, drop0trailing = TRUE, scientific = 3, ...
  ){
    format(
      pretty_round(
        x, round_digits, round_digits
      ), drop0trailing = drop0trailing,
      scientific = scientific, ...
    )
  }
  abbr <- function(x, min, left = FALSE){
    abbreviate(x, minlength = min, named = FALSE,
               method = if (left) "left.kept" else "both.sides")
  }
  format_num_cols <- function(data, fun){
    for (i in seq_along(data)){
      if (is.numeric(data[[i]])){
        data[[i]] <- fun(data[[i]])
      }
    }
    data
  }
  for (i in seq_along(temp)){
    # maybe_tbl <- package_has_function("tibble", "as_tibble")
    # if (maybe_tbl){
    #   as_tibble <- get("as_tibble", asNamespace("tibble"), inherits = FALSE)
    # }
    if (inherits(temp[[i]], "data.frame")){
      # if (maybe_tbl){
      #   temp[[i]] <- as_tibble(format_num_cols(temp[[i]], pretty_round))
      # } else {
      temp[[i]] <- format_num_cols(temp[[i]], pretty_num)
      temp[[i]][["class"]] <- NULL
      # }
    }
  }
  cat(paste("obs:", temp$obs, "\ncols:", temp$cols), "\n")
  if (nrow(temp$logical)){
    cat("\n----- Logical -----\n")
    names(temp$logical) <- abbr(names(temp$logical), 8)
    temp$logical$col <- abbr(temp$logical$col, 14)
    print(temp$logical)
  }
  if (nrow(temp$numeric)){
    cat("\n----- Numeric -----\n")
    names(temp$numeric) <- abbr(names(temp$numeric), 8)
    temp$numeric$col <- abbr(temp$numeric$col, 14)
    print(temp$numeric)
  }
  if (nrow(temp$date)){
    cat("\n----- Dates -----\n")
    names(temp$date) <- abbr(names(temp$date), 8)
    temp$date$col <- abbr(temp$date$col, 14)
    print(temp$date)
  }
  if (nrow(temp$datetime)){
    # An overview list contains a 'min' & 'max' variable of date-times
    # This is UTC because R can't handle a date-time with multiple time-zones
    # And so we want to print it in local-time
    datetime_chr_min <- character(nrow(temp$datetime))
    datetime_chr_max <- character(nrow(temp$datetime))
    mins <- temp[["datetime"]][["min"]]
    maxs <- temp[["datetime"]][["max"]]
    tzones <- temp[["datetime"]][["tzone"]]
    for (i in seq_len(nrow(temp$datetime))){
      datetime_chr_min[i] <- format(mins[i], tz = tzones[i])
      datetime_chr_max[i] <- format(maxs[i], tz = tzones[i])
    }
    temp$datetime$min <- datetime_chr_min
    temp$datetime$max <- datetime_chr_max
    cat("\n----- Date-Times -----\n")
    names(temp$datetime) <- abbr(names(temp$datetime), 8)
    temp$datetime$col <- abbr(temp$datetime$col, 14)
    print(temp$datetime)
  }
  if (nrow(temp$time_series)){
    cat("\n----- Time-Series -----\n")
    temp$time_series$growth_rate <- paste0(
      pretty_num(
        (x$time_series$growth_rate - 1) * 100
      ), "%"
    )
    names(temp$time_series) <- abbr(names(temp$time_series), 8)
    temp$time_series$col <- abbr(temp$time_series$col, 14)
    print(temp$time_series)
  }
  if (nrow(temp$categorical)){
    cat("\n----- Categorical -----\n")
    names(temp$categorical) <- abbr(names(temp$categorical), 8)
    temp$categorical$col <- abbr(temp$categorical$col, 14)
    print(temp$categorical)
  }
  if (nrow(temp$other)){
    cat("\n----- Other -----\n")
    names(temp$other) <- abbr(names(temp$other), 8)
    temp$other$col <- abbr(temp$other$col, 14)
    print(temp$other)
  }
  invisible(x)
}

# Helpers -----------------------------------------------------------------


n_unique <- function(x, na_rm = FALSE){
  out <- collapse::fnunique(x)
  if (na_rm){
    out <- out - any_na(x, recursive = FALSE)
  }
  out
}
prop_missing <- function(x, recursive = TRUE){
  if (recursive){
    N <- unlisted_length(x)
  } else {
    N <- cpp_vec_length(x)
  }
  na_count(x, recursive = recursive) / N
}
prop_complete <- function(x, recursive = TRUE){
  1 - prop_missing(x, recursive = recursive)
}
transform_all <- function(data, .fn, .cols = names(data)){
  out <- unclass(data)
  for (col in .cols){
    out[[col]] <- .fn(out[[col]])
  }
  class(out) <- class(data)
  out
}
summarise_all <- function(data, .fn, size = 1){
  out <- unclass(sset(data, seq_len(size)))
  for (col in names(out)){
    out[[col]] <- .fn(data[[col]])
  }
  attr(out, "row.names") <- .set_row_names(size)
  class(out) <- class(data)
  out
}
pluck_row <- function(data, i = 1){
  unlist(sset(data, i), recursive = FALSE)
}

# Taken from skimr::skim with modifications
spark_bar <- function(x){
  bars <- intToUtf8(c(9601L, 9602L, 9603L, 9605L, 9606L, 9607L),
                    multiple = TRUE)
  bar_codes <- bin(
    x, seq.int(0, to = 1, length.out = length(bars) + 1L),
    left_closed = TRUE, include_oob = TRUE,
    include_endpoint = TRUE
  )
  bar_codes[bar_codes == 0L] <- NA_integer_
  out <- bars[bar_codes]
  paste0(out, collapse = "")
}

inline_hist <- function(x, n_bins = 5L){
  if (length(x) < 1L) {
    return("")
  }
  if (is.infinite(max(abs(collapse::frange(x, na.rm = TRUE))))){
    x[which_(is.infinite(x))] <- NA
  }
  n_nas <- na_count(x)
  all_na <- n_nas == length(x)

  if (all_na){
    return("")
  }
  # All zeros?
  if (val_count(x, 0) == (length(x) - n_nas)){
    x <- x + 1
  }
  hist_dt <- tabulate(
    bin(x, r_cut_breaks(x, n_bins), left_closed = FALSE),
    nbins = n_bins
  )
  hist_dt <- hist_dt / max(hist_dt)
  spark_bar(hist_dt)
}

# Efficient summary statistics
numeric_summary <- function(x){
  n <- length(x)
  if (n == 0){
    return(
      fast_df(
        n_missing = NA_integer_,
        p_complete = NA_real_,
        n_unique = NA_integer_,
        mean = NA_real_,
        p0 = NA_real_,
        p25 = NA_real_,
        p50 = NA_real_,
        p75 = NA_real_,
        p100 = NA_real_,
        iqr = NA_real_,
        sd = NA_real_
      )
    )
  }
  n_missing <- na_count(x)
  n_unique <- n_unique(x) - (n_missing > 0L)
  p_complete <- 1 - (n_missing / n)
  sum <- sum(x, na.rm = TRUE)
  mean <- sum / (n - n_missing)

  if (!is.integer(x)){
    x <- as.double(x)
  }
  # o <- order(x, na.last = TRUE)
  # percentiles <- c(
  #   x[o[1]], # min
  #   collapse::fnth(x, 0.25, na.rm = TRUE, o = o),
  #   collapse::fnth(x, 0.5, na.rm = TRUE, o = o),
  #   collapse::fnth(x, 0.75, na.rm = TRUE, o = o),
  #   x[o[length(o) - n_missing]] # max
  # )
  rng <- collapse::frange(x, na.rm = TRUE)
  percentiles <- c(
    rng[1],
    collapse::fnth(x, 0.25, na.rm = TRUE),
    collapse::fnth(x, 0.5, na.rm = TRUE),
    collapse::fnth(x, 0.75, na.rm = TRUE),
    rng[2]
  )
  # percentiles <- collapse::fquantile(
  #   x, probs = c(0, 0.25, 0.5, 0.75, 1),
  #   na.rm = TRUE
  # )
  var <- var_sum_squared_diff(x, mean) / (n - n_missing - 1)
  p0 <- percentiles[1]
  p25 <- percentiles[2]
  p50 <- percentiles[3]
  p75 <- percentiles[4]
  p100 <- percentiles[5]
  sd <- sqrt(var)
  iqr <- p75 - p25

  fast_df(
    n_missing = n_missing,
    p_complete = p_complete,
    n_unique = n_unique,
    mean = mean,
    p0 = p0,
    p25 = p25,
    p50 = p50,
    p75 = p75,
    p100 = p100,
    iqr = iqr,
    sd = sd
  )
}

# A custom rounding that uses decimal places for large numbers
# and significant digits for small numbers
pretty_round <- function(x, decimal_digits = 2, signif_digits = 2, threshold = 1){
  if (!is.numeric(x)){
    stop("x must be numeric")
  }
  if (!is.integer(x) && is.numeric(x)){
    x <- cheapr_if_else(
      abs(x) < threshold,
      signif(x, signif_digits),
      round(x, decimal_digits)
    )
  }
  x
}
# ts_growth_rate_per_period <- function(x){
#   time <- as.vector(stats::time(x))
#   freq <- stats::frequency(x)
#   a <- utils::head(x, n = 1)
#   b <- utils::tail(x, n = 1)
#   cpp_growth_rate(as.double(a), as.double(b), length(time) / freq)
# }

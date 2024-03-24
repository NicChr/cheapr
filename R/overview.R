#' An alternative to `summary()` inspired by the skimr package
#'
#' @description
#' A cheaper `summary()` function, designed for larger data.
#'
#' @param x A vector or data frame.
#' @param hist Should in-line histograms be returned? Default is `FALSE`.
#' @param digits How many decimal places should the summary statistics be
#' printed as? Default is 2.
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
#' @rdname overview
#' @export
overview <- function(x, hist = FALSE, digits = getOption("cheapr.digits", 2)){
 UseMethod("overview")
}
#' @rdname overview
#' @export
overview.default <- function(x, hist = FALSE, digits = getOption("cheapr.digits", 2)){
  options(cheapr.digits = digits)
  overview(list_as_df(list(x = x)), hist = hist)
}
#' @rdname overview
#' @export
overview.logical <- function(x, hist = FALSE, digits = getOption("cheapr.digits", 2)){
  options(cheapr.digits = digits)
  overview(list_as_df(list(x = as.logical(x))), hist = hist)
}
#' @rdname overview
#' @export
overview.numeric <- function(x, hist = FALSE, digits = getOption("cheapr.digits", 2)){
  options(cheapr.digits = digits)
  out <- overview(list_as_df(list(x = as.numeric(x))), hist = hist)
  out$cols <- NA_integer_
  out
}
#' @rdname overview
#' @export
overview.character <- function(x, hist = FALSE, digits = getOption("cheapr.digits", 2)){
  options(cheapr.digits = digits)
  out <- overview(list_as_df(list(x = as.character(x))), hist = hist)
  out$cols <- NA_integer_
  out
}
#' @rdname overview
#' @export
overview.factor <- function(x, hist = FALSE, digits = getOption("cheapr.digits", 2)){
  options(cheapr.digits = digits)
  out <- overview(list_as_df(list(x = as.factor(x))), hist = hist)
  out$cols <- NA_integer_
  out
}
#' @rdname overview
#' @export
overview.Date <- function(x, hist = FALSE, digits = getOption("cheapr.digits", 2)){
  options(cheapr.digits = digits)
  out <- overview(list_as_df(list(x = as.Date(x))), hist = hist)
  out$cols <- NA_integer_
  out
}
#' @rdname overview
#' @export
overview.POSIXt <- function(x, hist = FALSE, digits = getOption("cheapr.digits", 2)){
  options(cheapr.digits = digits)
  out <- overview(list_as_df(list(x = as.POSIXct(x))), hist = hist)
  out$cols <- NA_integer_
  out
}
#' @rdname overview
#' @export
overview.ts <- function(x, hist = FALSE, digits = getOption("cheapr.digits", 2)){
  options(cheapr.digits = digits)
  out <- overview(transform_all(as.data.frame(x), as.numeric), hist = hist)
  out$time_series <- out$numeric
  out$numeric <- sset(out$numeric, 0)
  out$time_series$class <- class(x)[1]
  out
}
#' @rdname overview
#' @export
overview.zoo <- overview.ts
#' @rdname overview
#' @export
overview.data.frame <- function(x, hist = FALSE, digits = getOption("cheapr.digits", 2)){
  options(cheapr.digits = digits)
  check_is_df(x)
  N <- nrow(x)
  num_cols <- ncol(x)
  skim_df <- x
  data_nms <- names(skim_df)
  col_classes <- vapply(skim_df, function(x) utils::tail(class(x), n = 1), "")
  out <- list_as_df(enframe_(col_classes, name = "col", value = "class"))
  chr_vars <- data_nms[vapply(skim_df, is.character, FALSE,
                              USE.NAMES = FALSE)]
  if (length(chr_vars) > 0L) {
    for (chr in chr_vars){
      skim_df[[chr]] <- factor_(skim_df[[chr]], ordered = TRUE)
    }
  }
  lgl_vars <- data_nms[vapply(skim_df, is.logical, FALSE)]
  num_vars <- data_nms[vapply(skim_df, function(x) inherits(x,
                                                            c("integer", "numeric")), FALSE)]
  date_vars <- data_nms[vapply(skim_df, function(x) inherits(x, "Date"), FALSE)]
  datetime_vars <- data_nms[vapply(skim_df, function(x) inherits(x, "POSIXt"),  FALSE)]
  ts_vars <- data_nms[vapply(skim_df, function(x) inherits(x, c("ts", "zoo")),  FALSE)]
  cat_vars <- data_nms[vapply(skim_df, is.factor, FALSE)]
  other_vars <- setdiff_(data_nms, c(lgl_vars, num_vars, date_vars,
                                     datetime_vars, ts_vars, cat_vars))

# Logical -----------------------------------------------------------------

  lgl_data <- df_select(skim_df, lgl_vars)
  which_lgl <- which_in(out[["col"]], lgl_vars)
  lgl_out <- sset(out, which_lgl)
  value_size <- min(length(which_lgl), 1L)
  lgl_out <- df_add_cols(lgl_out, list(n_missing = NA_integer_[value_size]))
  lgl_out <- df_add_cols(lgl_out, list(p_complete = NA_real_[value_size]))
  lgl_out <- df_add_cols(lgl_out, list(n_true = NA_integer_[value_size],
                                       n_false = NA_integer_[value_size]))
  lgl_out <- df_add_cols(lgl_out, list(p_true = NA_real_[value_size]))
  if (N > 0L && length(which_lgl) > 0) {
    lgl_out$n_missing <- pluck_row(summarise_all(lgl_data, num_na), 1)
    lgl_out$p_complete <- pluck_row(summarise_all(lgl_data, prop_complete), 1)
    lgl_out$n_true <- pluck_row(summarise_all(lgl_data, function(x) sum(x, na.rm = TRUE)), 1)
    lgl_out$n_false <- N - lgl_out[["n_missing"]] - lgl_out[["n_true"]]
    lgl_out$p_true <- lgl_out[["n_true"]] / (N - lgl_out[["n_missing"]])
  }

# Numeric -----------------------------------------------------------------

  num_data <- df_select(skim_df, num_vars)
  which_num <- which_in(out[["col"]], num_vars)
  num_out <- sset(out, which_num)
  value_size <- min(length(which_num), 1L)
  num_out <- df_add_cols(num_out, list(n_missing = NA_integer_[value_size]))
  num_out <- df_add_cols(num_out, list(p_complete = NA_real_[value_size]))
  num_out <- df_add_cols(num_out, list(n_unique = NA_integer_[value_size]))
  num_out <- df_add_cols(num_out, stats::setNames(
    new_list(8, default = NA_real_[value_size]),
    c("mean", "p0", "p25", "p50", "p75", "p100", "iqr", "sd")
  ))
  if (hist){
    num_out$hist <- NA_character_[value_size]
  }
  if (N > 0L && length(which_num) > 0) {
    num_out$n_missing <- pluck_row(summarise_all(num_data, num_na), 1)
    num_out$p_complete <- pluck_row(summarise_all(num_data, prop_complete), 1)
    num_out$n_unique <- pluck_row(summarise_all(num_data, n_unique), 1)
    num_out$n_unique <- num_out$n_unique - (num_out$n_missing > 0L)
    num_out$mean <- pluck_row(summarise_all(num_data, collapse::fmean), 1)
    num_out$p0 <- pluck_row(summarise_all(num_data, collapse::fmin), 1)
    num_out$p25 <- pluck_row(summarise_all(
      num_data, function(x) collapse::fnth(x, 0.25)
    ), 1)
    num_out$p50 <- pluck_row(summarise_all(
      num_data, function(x) collapse::fnth(x, 0.5)
    ), 1)
    num_out$p75 <- pluck_row(summarise_all(
      num_data, function(x) collapse::fnth(x, 0.75)
    ), 1)
    num_out$p100 <- pluck_row(summarise_all(num_data, collapse::fmax), 1)
    num_out$iqr <- num_out$p75 - num_out$p25
    num_out$sd <- pluck_row(summarise_all(num_data, collapse::fsd), 1)
    if (hist){
      num_out$hist <- pluck_row(summarise_all(
        num_data, inline_hist
      ), 1)
    }
  }

# Dates -------------------------------------------------------------------

  date_data <- df_select(skim_df, date_vars)
  which_date <- which_in(out[["col"]], date_vars)
  date_out <- sset(out, which_date)
  value_size <- min(length(which_date), 1L)
  date_out <- df_add_cols(date_out, list(n_missing = NA_integer_[value_size]))
  date_out <- df_add_cols(date_out, list(p_complete = NA_real_[value_size]))
  date_out <- df_add_cols(date_out, list(n_unique = NA_integer_[value_size]))
  date_out <- df_add_cols(date_out, list(min = .Date(NA_real_[value_size]),
                                         max = .Date(NA_real_[value_size])))
  if (N > 0L && length(which_date) > 0) {
    date_out$n_missing <- pluck_row(summarise_all(date_data, num_na), 1)
    date_out$p_complete <- pluck_row(summarise_all(date_data, prop_complete), 1)
    date_out$n_unique <- pluck_row(summarise_all(date_data, n_unique), 1)
    date_out$n_unique <- date_out$n_unique - (date_out$n_missing > 0L)
    date_out$min <- pluck_row(summarise_all(date_data, collapse::fmin), 1)
    date_out$max <- pluck_row(summarise_all(date_data, collapse::fmax), 1)
    class(date_out$min) <- "Date"
    class(date_out$max) <- "Date"
  }


# Date-Times --------------------------------------------------------------

  datetime_data <- df_select(skim_df, datetime_vars)
  datetime_data <- transform_all(datetime_data, as.POSIXct)
  which_datetime <- which_in(out[["col"]], datetime_vars)
  datetime_out <- sset(out, which_datetime)
  value_size <- min(length(which_datetime), 1L)
  datetime_out <- df_add_cols(datetime_out, list(n_missing = NA_integer_[value_size]))
  datetime_out <- df_add_cols(datetime_out, list(p_complete = NA_real_[value_size]))
  datetime_out <- df_add_cols(datetime_out, list(n_unique = NA_integer_[value_size]))
  datetime_out <- df_add_cols(datetime_out, list(tzone = NA_character_[value_size]))
  datetime_out <- df_add_cols(datetime_out, list(min = .POSIXct(NA_real_[value_size]),
                                                 max = .POSIXct(NA_real_[value_size])))
  if (N > 0L && length(which_datetime) > 0) {
    datetime_out$n_missing <- pluck_row(summarise_all(datetime_data, num_na), 1)
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

  ts_data <- df_select(skim_df, ts_vars)
  which_ts <- which_in(out[["col"]], ts_vars)
  ts_out <- sset(out, which_ts)
  if (N > 0L && length(which_ts) > 0) {
    ts_overviews <- new_list(nrow(ts_out))
    for (i in seq_along(ts_overviews)){
      ts_overviews[[i]] <- overview(ts_data[[ts_out[["col"]][i]]], hist = hist)$time_series
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

  cat_data <- df_select(skim_df, cat_vars)
  which_cat <- which_in(out[["col"]], cat_vars)
  cat_out <- sset(out, which_cat)
  value_size <- min(length(which_cat), 1L)
  cat_out <- df_add_cols(cat_out, list(n_missing = NA_integer_[value_size]))
  cat_out <- df_add_cols(cat_out, list(p_complete = NA_real_[value_size]))
  cat_out <- df_add_cols(cat_out, list(n_unique = NA_integer_[value_size]))
  cat_out <- df_add_cols(cat_out, list(n_levels = NA_integer_[value_size]))
  cat_out <- df_add_cols(cat_out, list(min = NA_character_[value_size],
                                       max = NA_character_[value_size]))
  if (N > 0L && length(which_cat) > 0) {
    cat_out$n_missing <- pluck_row(summarise_all(cat_data, num_na), 1)
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

  other_data <- df_select(skim_df, other_vars)
  which_other <- which_in(out[["col"]], other_vars)
  other_out <- sset(out, which_other)
  value_size <- min(length(which_other), 1L)
  other_out <- df_add_cols(other_out, list(n_missing = NA_integer_[value_size]))
  other_out <- df_add_cols(other_out, list(p_complete = NA_real_[value_size]))
  other_out <- df_add_cols(other_out, list(n_unique = NA_integer_[value_size]))
  if (N > 0L && length(which_other) > 0) {
    other_out$n_missing <- pluck_row(summarise_all(
      other_data, function(x) num_na(x, recursive = FALSE)
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
    logical = lgl_out,
    numeric = num_out,
    date = date_out,
    datetime = datetime_out,
    time_series = ts_out,
    categorical = cat_out,
    other = other_out
  )
  class(out) <- c("overview", "list")
  out
}
#' @export
print.overview <- function(x, max = NULL, digits = getOption("cheapr.digits", 2), ...){
  # max_rows <- getOption("tibble.print_max", 20)
  # max_cols <- getOption("tibble.width", NULL)
  # max_extra_cols <- getOption("tibble.max_extra_cols", 100)
  # options(tibble.print_max = 10)
  # options(tibble.width = 100)
  # options(tibble.max_extra_cols = 10)
  cat(paste("obs:", x$obs, "\ncols:", x$cols), "\n")
  # for (data_type in names(x)[-(1:2)]){
  #   if (nrow(x[[data_type]])){
  #     cat(paste("\n-----", data_type, "-----\n"))
  #     print(x[[data_type]])
  #   }
  # }
  if (nrow(x$logical)){
    x$logical$p_complete <- pretty_num(round(x$logical$p_complete, digits))
    cat("\n----- Logical -----\n")
    print(x$logical)
  }
  if (nrow(x$numeric)){
    x$numeric$p_complete <- pretty_num(round(x$numeric$p_complete, digits))
    x$numeric$mean <- pretty_num(round(x$numeric$mean, digits))
    x$numeric$p0 <- pretty_num(round(x$numeric$p0, digits))
    x$numeric$p25 <- pretty_num(round(x$numeric$p25, digits))
    x$numeric$p50 <- pretty_num(round(x$numeric$p50, digits))
    x$numeric$p75 <- pretty_num(round(x$numeric$p75, digits))
    x$numeric$p100 <- pretty_num(round(x$numeric$p100, digits))
    x$numeric$iqr <- pretty_num(round(x$numeric$iqr, digits))
    x$numeric$sd <- pretty_num(round(x$numeric$sd, digits))
    cat("\n----- Numeric -----\n")
    print(x$numeric)
  }
  if (nrow(x$date)){
    x$date$p_complete <- pretty_num(round(x$date$p_complete, digits))
    cat("\n----- Dates -----\n")
    print(x$date)
  }
  if (nrow(x$datetime)){
    x$datetime$p_complete <- pretty_num(round(x$datetime$p_complete, digits))
    # An overview list contains a 'min' & 'max' variable of date-times
    # This is UTC because R can't handle a date-time with multiple time-zones
    # And so we want to print it in local-time
    datetime_chr_min <- character(nrow(x$datetime))
    datetime_chr_max <- character(nrow(x$datetime))
    mins <- x[["datetime"]][["min"]]
    maxs <- x[["datetime"]][["max"]]
    tzones <- x[["datetime"]][["tzone"]]
    for (i in seq_len(nrow(x$datetime))){
      datetime_chr_min[i] <- format(mins[i], tz = tzones[i])
      datetime_chr_max[i] <- format(maxs[i], tz = tzones[i])
    }
    x$datetime$min <- datetime_chr_min
    x$datetime$max <- datetime_chr_max
    cat("\n----- Date-Times -----\n")
    print(x$datetime)
  }
  if (nrow(x$time_series)){
    x$time_series$p_complete <- pretty_num(round(x$time_series$p_complete, digits))
    x$time_series$mean <- pretty_num(round(x$time_series$mean, digits))
    x$time_series$p0 <- pretty_num(round(x$time_series$p0, digits))
    x$time_series$p25 <- pretty_num(round(x$time_series$p25, digits))
    x$time_series$p50 <- pretty_num(round(x$time_series$p50, digits))
    x$time_series$p75 <- pretty_num(round(x$time_series$p75, digits))
    x$time_series$p100 <- pretty_num(round(x$time_series$p100, digits))
    x$time_series$iqr <- pretty_num(round(x$time_series$iqr, digits))
    x$time_series$sd <- pretty_num(round(x$time_series$sd, digits))
    cat("\n----- Time-Series -----\n")
    print(x$time_series)
  }
  if (nrow(x$categorical)){
    x$categorical$p_complete <- pretty_num(round(x$categorical$p_complete, digits))
    cat("\n----- Categorical -----\n")
    print(x$categorical)
  }
  if (nrow(x$other)){
    x$other$p_complete <- pretty_num(round(x$other$p_complete, digits))
    cat("\n----- Other -----\n")
    print(x$other)
  }
  # options(tibble.print_max = max_rows)
  # options(tibble.width = max_cols)
  # options(tibble.max_extra_cols = max_extra_cols)
  invisible(x)
}
### Helpers

n_unique <- function(x, na_rm = FALSE){
  out <- collapse::fnunique(x)
  if (na_rm){
    out <- out - any_na(x, recursive = FALSE)
  }
  out
}
prop_complete <- function(x, recursive = TRUE){
  if (recursive){
    N <- unlisted_length(x)
  } else {
    N <- length(x)
  }
  1 - (num_na(x, recursive = recursive) / N)
}
transform_all <- function(data, .fn){
  for (col in names(data)){
    data[[col]] <- .fn(data[[col]])
  }
  data
}
summarise_all <- function(data, .fn, size = 1){
  out <- sset(data, seq_len(size))
  attr(out, "row.names") <- .set_row_names(size)
  for (col in names(out)){
    out[[col]] <- .fn(data[[col]])
  }
  out
}
pluck_row <- function(data, i = 1){
  unlist(data[i, ], recursive = FALSE)
}

# Taken from skimr::skim with modifications
spark_bar <- function(x){
  bars <- intToUtf8(c(9601L, 9602L, 9603L, 9605L, 9606L, 9607L),
                    multiple = TRUE)
  bar_codes <- findInterval(
    x, vec = seq.int(0, to = 1, length.out = length(bars) + 1L),
    rightmost.closed = TRUE,
    left.open = FALSE, all.inside = FALSE
  )
  bar_codes[bar_codes == 0L] <- NA_integer_
  out <- bars[bar_codes]
  paste0(out, collapse = "")
}

inline_hist <- function(x, n_bins = 5L){
  if (length(x) < 1L) {
    return(" ")
  }
  if (is.infinite(max(abs(collapse::frange(x, na.rm = TRUE))))) {
    x[is.infinite(x)] <- NA
  }
  if (all_na(x)) {
    return(" ")
  }
  if (allv2(na_rm(x), 0)) {
    x <- x + 1
  }
  hist_dt <- tabulate(cut_numeric(x, n_bins, labels = FALSE),
                      nbins = n_bins)
  hist_dt <- hist_dt / max(hist_dt)
  spark_bar(hist_dt)
}
pretty_num <- function(x, scientific = FALSE, ...){
  prettyNum(x, scientific = scientific, ...)
}

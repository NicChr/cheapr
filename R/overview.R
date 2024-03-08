#' An alternative to `summary()` inspired by the skimr package
#'
#' @description
#' A cheaper `summary()` function, designed for larger data.
#'
#' @param x A vector or data frame.
#' @param hist Should in-line histograms be returned? Default is `FALSE`.
#'
#' @returns
#' `overview(x)` returns a 1-row data frame unless
#' `x` is a data frame, in which a list of data frames is returned.
#' Key summary statistics are reported in each data frame.
#'
#' @export
overview <- function(x, hist = FALSE){
 UseMethod("overview")
}
#' @export
overview.default <- function(x, hist = FALSE){
  warning("Unsure how to calculate summary for x, falling back to character")
  out <- as.character(x)
  out <- overview(list_as_df(list(x = out)), hist = hist)$categorical
  names(out)[1] <- "length"
  out[[1]] <- length(x)
  out[[2]] <- utils::tail(class(x), n = 1)
  out
}
#' @export
overview.logical <- function(x, hist = FALSE){
  out <- overview(list_as_df(list(x = x)), hist = hist)$logical
  names(out)[1] <- "length"
  out[[1]] <- length(x)
  out
}
#' @export
overview.numeric <- function(x, hist = FALSE){
  out <- overview(list_as_df(list(x = x)), hist = hist)$numeric
  names(out)[1] <- "length"
  out[[1]] <- length(x)
  out
}
#' @export
overview.character <- function(x, hist = FALSE){
  out <- overview(list_as_df(list(x = x)), hist = hist)$categorical
  names(out)[1] <- "length"
  out[[1]] <- length(x)
  out
}
#' @export
overview.Date <- function(x, hist = FALSE){
  out <- overview(list_as_df(list(x = x)), hist = hist)$date
  names(out)[1] <- "length"
  out[[1]] <- length(x)
  out
}
#' @export
overview.POSIXt <- function(x, hist = FALSE){
  out <- overview(list_as_df(list(x = as.POSIXct(x))), hist = hist)$datetime
  names(out)[1] <- "length"
  out[[1]] <- length(x)
  out[[2]] <- utils::tail(class(x), n = 1)
  out
}
#' @export
overview.data.frame <- function(x, hist = FALSE){
  check_is_df(x)
  N <- nrow(x)
  num_cols <- ncol(x)
  skim_df <- df_as_tbl(list_as_df(x))
  data_nms <- names(skim_df)
  col_classes <- vapply(skim_df, function(x) utils::tail(class(x), n = 1), "")
  out <- df_as_tbl(enframe_(col_classes, name = "col", value = "class"))
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
  cat_vars <- data_nms[vapply(skim_df, is.factor, FALSE)]
  other_vars <- setdiff(data_nms, c(lgl_vars, num_vars, date_vars,
                                    datetime_vars, cat_vars))
  if (length(other_vars) > 0) {
    warning(paste0("Unsure how to calculate summaries for these variables: \n",
                   paste(other_vars, collapse = "\n")), "\n\nFalling back to character")
    for (var in other_vars) {
      skim_df[[var]] <- factor_(as.character(skim_df[[var]]), ordered = TRUE)
    }
    cat_vars <- c(cat_vars, other_vars)
    cat_vars <- data_nms[sort(match(cat_vars, data_nms))]
  }

# Logical -----------------------------------------------------------------

  lgl_data <- df_select(skim_df, lgl_vars)
  which_lgl <- which_in(out[["col"]], lgl_vars)
  lgl_out <- out[which_lgl, , drop = FALSE]
  lgl_out <- df_add_cols(lgl_out, list(n_missing = NA_integer_))
  lgl_out <- df_add_cols(lgl_out, list(p_complete = NA_real_))
  lgl_out <- df_add_cols(lgl_out, list(n_true = NA_integer_, n_false = NA_integer_))
  lgl_out <- df_add_cols(lgl_out, list(p_true = NA_real_))
  if (N > 0L && length(which_lgl) > 0) {
    lgl_out$n_missing <- pluck_row(summarise_all(lgl_data, num_na), 1)
    lgl_out$p_complete <- pluck_row(summarise_all(lgl_data, prop_complete), 1)
    lgl_out$n_true <- pluck_row(summarise_all(lgl_data, function(x) sum(x, na.rm = TRUE)), 1)
    lgl_out$n_false <- N - lgl_out[["n_missing"]] - lgl_out[["n_true"]]
    lgl_out$p_true <- lgl_out[["n_true"]] / (N - lgl_out[["n_missing"]])
  }
  lgl_out <- df_as_tbl(lgl_out)

# Numeric -----------------------------------------------------------------

  num_data <- df_select(skim_df, num_vars)
  which_num <- which_in(out[["col"]], num_vars)
  num_out <- out[which_num, , drop = FALSE]
  num_out <- df_add_cols(num_out, list(n_missing = NA_integer_))
  num_out <- df_add_cols(num_out, list(p_complete = NA_real_))
  num_out <- df_add_cols(num_out, list(n_unique = NA_integer_))
  num_out <- df_add_cols(num_out, stats::setNames(as.list(rep_len(NA_real_,
                                                  8)), c("mean", "p0", "p25", "p50", "p75", "p100", "iqr",
                                                         "sd")))
  if (hist){
    num_out$hist <- NA_character_
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
  num_out <- df_as_tbl(num_out)

# Dates -------------------------------------------------------------------

  date_data <- df_select(skim_df, date_vars)
  which_date <- which_in(out[["col"]], date_vars)
  date_out <- out[which_date, , drop = FALSE]
  date_out <- df_add_cols(date_out, list(n_missing = NA_integer_))
  date_out <- df_add_cols(date_out, list(p_complete = NA_real_))
  date_out <- df_add_cols(date_out, list(n_unique = NA_integer_))
  date_out <- df_add_cols(date_out, list(min = .Date(NA_real_),
                                         max = .Date(NA_real_)))
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
  date_out <- df_as_tbl(date_out)


# Date-Times --------------------------------------------------------------

  datetime_data <- df_select(skim_df, datetime_vars)
  datetime_data <- transform_all(datetime_data, as.POSIXct)
  which_datetime <- which_in(out[["col"]], datetime_vars)
  datetime_out <- out[which_datetime, , drop = FALSE]
  datetime_out <- df_add_cols(datetime_out, list(n_missing = NA_integer_))
  datetime_out <- df_add_cols(datetime_out, list(p_complete = NA_real_))
  datetime_out <- df_add_cols(datetime_out, list(n_unique = NA_integer_))
  datetime_out <- df_add_cols(datetime_out, list(min = .POSIXct(NA_real_),
                                                 max = .POSIXct(NA_real_)))
  if (N > 0L && length(which_datetime) > 0) {
    datetime_out$n_missing <- pluck_row(summarise_all(datetime_data, num_na), 1)
    datetime_out$p_complete <- pluck_row(summarise_all(datetime_data, prop_complete), 1)
    datetime_out$n_unique <- pluck_row(summarise_all(datetime_data, n_unique), 1)
    datetime_out$n_unique <- datetime_out$n_unique - (datetime_out$n_missing > 0L)
    datetime_out$min <- pluck_row(summarise_all(datetime_data, collapse::fmin), 1)
    datetime_out$max <- pluck_row(summarise_all(datetime_data, collapse::fmax), 1)
    datetime_out$min <- .POSIXct(datetime_out$min, tz = "UTC")
    datetime_out$max <- .POSIXct(datetime_out$max, tz = "UTC")
  }
  datetime_out <- df_as_tbl(datetime_out)


# Categorical -------------------------------------------------------------

  cat_data <- df_select(skim_df, cat_vars)
  which_cat <- which_in(out[["col"]], cat_vars)
  cat_out <- out[which_cat, , drop = FALSE]
  cat_out <- df_add_cols(cat_out, list(n_missing = NA_integer_))
  cat_out <- df_add_cols(cat_out, list(p_complete = NA_real_))
  cat_out <- df_add_cols(cat_out, list(n_unique = NA_integer_))
  cat_out <- df_add_cols(cat_out, stats::setNames(as.list(rep_len(NA_character_, 2)),
                                            c("min", "max")))
  if (N > 0L && length(which_cat) > 0) {
    cat_out$n_missing <- pluck_row(summarise_all(cat_data, num_na), 1)
    cat_out$p_complete <- pluck_row(summarise_all(cat_data, prop_complete), 1)
    cat_out$n_unique <- pluck_row(summarise_all(cat_data, n_unique), 1)
    cat_out$n_unique <- cat_out$n_unique - (cat_out$n_missing > 0L)
    cat_out$min <- pluck_row(summarise_all(cat_data, collapse::fmin), 1)
    cat_out$max <- pluck_row(summarise_all(cat_data, collapse::fmax), 1)
  }
  cat_out <- df_as_tbl(cat_out)
  list(
    nrow = N, ncol = num_cols,
    logical = lgl_out,
    numeric = num_out,
    date = date_out,
    datetime = datetime_out,
    categorical = cat_out
  )
}

### Helpers
n_unique <- function(x){
  collapse::fnunique(x)
}
prop_complete <- function(x){
  N <- unlisted_length(x)
  1 - (num_na(x) / N)
}
transform_all <- function(data, .fn){
  for (col in names(data)){
    data[[col]] <- .fn(data[[col]])
  }
  data
}
summarise_all <- function(data, .fn, size = 1){
  out <- data[seq_len(size), , drop = FALSE]
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

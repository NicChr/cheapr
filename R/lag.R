#' Lagged operations.
#'
#' @description
#' Fast lags and leads
#'
#' @param x A vector or data frame.
#' @param n Number of lags. Negative values are accepted.
#' @param fill Value used to fill first n values. Default is `NA`.
#' @param set Should x be updated by reference? If `TRUE` no copy is made and
#' x is updated in place. The default is `FALSE`.
#' @param recursive Should list elements be lagged as well?
#' If `TRUE`, this is useful for data frames and will return row lags.
#' If `FALSE` this will return a plain lagged list.
#'
#' @returns
#' A lagged object the same size as x.
#'
#' @examples
#' library(cheapr)
#' library(bench)
#'
#' # A use-case for data.table
#'
#' df <- data.frame(x = 1:10^5)
#'
#' # Lag these behind by 3 rows
#' sset(lag_(df, 3, set = TRUE), 1:10)
#'
#' df$x[1:10] # x variable was updated by reference!
#'
#' # The above can be used naturally in data.table to lag data
#' # without any copies
#'
#' # To perform regular R row lags, just make sure set is `FALSE`
#'
#' sset(lag_(as.data.frame(EuStockMarkets), 5), 1:10)
#' @export
lag_ <- function(x, n = 1, fill = NULL, set = FALSE, recursive = TRUE){
  .Call(`_cheapr_cpp_lag`, x, n, fill, set, recursive)
}

#' Lagged operations.
#'
#' @description
#' Fast lags and leads optionally using dynamic vectorised lags, ordering and
#' run lengths.
#'
#' @param x A vector or data frame.
#' @param n Number of lags. Negative values are accepted. \cr
#' `lag2_` accepts a vector of dynamic lags and leads
#' which gets recycled to the length of x.
#' @param fill Value used to fill first n values. Default is `NA`.
#' @param set Should x be updated by reference? If `TRUE` no copy is made and
#' x is updated in place. The default is `FALSE`.
#' @param recursive Should list elements be lagged as well?
#' If `TRUE`, this is useful for data frames and will return row lags.
#' If `FALSE` this will return a plain lagged list.
#' @param order Optionally specify an ordering with which to apply the lags.
#' This is useful for example when applying lags chronologically using an unsorted
#' time variable.
#' @param run_lengths Optional integer vector of run lengths that defines the size of each
#' lag run. For example, supplying `c(5, 5)` applies lags to the first 5 elements and
#' then essentially resets the bounds and applies lags to the next 5 elements as if
#' they were an entirely separate and standalone vector. \cr
#' This is particularly useful in conjunction with the `order` argument
#' to perform a by-group lag. See the examples for details.
#'
#' @returns
#' A lagged object the same size as x.
#'
#' @details
#' For most applications, it is more efficient and recommended to use `lag_()`.
#' For anything that requires dynamic lags, lag by order of another variable,
#' or by-group lags, one can use `lag2_()`.
#'
#' ### `lag2_`
#'
#' `lag2_` is a generalised form of `lag_` that by default performs
#' simple lags and leads. \cr
#' It has 3 additional features but does not support updating by reference or
#' long vectors. \cr
#'
#' These extra features include:
#' * `n` - This shares the same name as the `n` argument in `lag_`
#' for consistency. The difference is that `lag_` accepts a
#' lag vector of length 1 whereas
#' this accepts a vector of dynamic lags allowing for
#' flexible combinations of variable sized lags and leads.
#' These are recycled to the length of the data and will always align
#' with the data, meaning that if you supply a custom `order` argument, this
#' ordering is applied both to `x` and the recycled lag vector `n` simultaneously.
#' * `order` - Apply lags in any order you wish. This can be useful for
#' reverse order lags, lags against unsorted time variables, and by-group lags.
#' * `run_lengths` - Specify the size of individual lag runs. For example, if you
#' specify `run_lengths = c(3, 4, 2)`, this will apply your lags to the first 3 elements and
#' then reset, applying lags to the next 4 elements, to reset again and apply lags to the
#' final 2 elements. Each time the reset occurs, it treats each run length sized 'chunk' as
#' a unique and separate vector. See the examples for a showcase.
#'
#' ### Table of differences between `lag_` and `lag2_`
#'
#' | Description | `lag_` | `lag2_` |
#' | :----: | :----: | :----: |
#' | Lags | Yes | Yes |
#' | Leads | Yes | Yes |
#' | Long vector support | Yes | No |
#' | Lag by reference | Yes | No |
#' | Dynamic vectorized lags | No | Yes |
#' | Data frame row lags | Yes | Yes |
#' | Alternative order lags | No | Yes |
#'
#' @examples
#' library(cheapr)
#' library(bench)
#'
#' # A use-case for data.table
#' # Adding 0 because can't update ALTREP by reference
#' df <- data.frame(x = 1:10^5 + 0L)
#'
#' # Normal data frame lag
#' sset(lag_(df), 1:10)
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
#'
#' # lag2_ is a generalised version of lag_ that allows
#' # for much more complex lags
#'
#' x <- 1:10
#'
#' # lag every 2nd element
#' lag2_(x, n = c(1, 0)) # lag vector is recycled
#'
#' # Explicit lags along x
#' lags <- lag_sequence(length(x), 1, partial = FALSE)
#' lag2_(x, n = lags)
#'
#' # Alternating lags and leads
#' lag2_(x, c(1, -1))
#'
#' # Lag only the 3rd element
#' lags <- integer(length(x))
#' lags[3] <- 1L
#' lag2_(x, lags)
#'
#' # lag in descending order (same as a lead)
#'
#' lag2_(x, order = 10:1)
#'
#' # lag that resets after index 5
#' lag2_(x, run_lengths = c(5, 5))
#'
#' # lag with a time index
#' years <- sample(2011:2020)
#' lag2_(x, order = order(years))
#'
#' # Example of how to do a cyclical lag
#' n <- length(x)
#' k <- min(3, n)
#' if (k >= 0){
#'   lag2_(x, c(rep(-n + k, k), rep(k, n - k)))
#' } else {
#'   lag2_(x, c(rep(k, n + k), rep(n + k, -k)))
#' }
#'
#' # As it turns out, we can do a grouped lag
#' # by supplying group sizes as run lengths and group order as the order
#'
#' set.seed(45)
#' g <- sample(c("a", "b"), 10, TRUE)
#'
#' # NOTE: collapse::flag will not work unless g is already sorted!
#' # This is not an issue with lag2_()
#' collapse::flag(x, g = g)
#' lag2_(x, order = order(g), run_lengths = collapse::GRP(g)$group.sizes)
#'
#' # For production code, we can of course make
#' # this more optimised by using collapse::radixorderv()
#' # Which calculates the order and group sizes all at once
#'
#' o <- collapse::radixorderv(g, group.sizes = TRUE)
#' lag2_(x, order = o, run_lengths = attr(o, "group.sizes"))
#'
#' # Let's finally wrap this up in a nice grouped-lag function
#'
#' grouped_lag <- function(x, n = 1, g = integer(length(x))){
#'   o <- collapse::radixorderv(g, group.sizes = TRUE, sort = FALSE)
#'   lag2_(x, n, order = o, run_lengths = attr(o, "group.sizes"))
#' }
#'
#' # And voila!
#' grouped_lag(x, g = g)
#'
#' # A method to extract this information from dplyr
#'
#' ## We can actually get this information easily from a `grouped_df` object
#' ## Uncomment the below code to run the implementation
#' # library(dplyr)
#' # library(timeplyr)
#' # eu_stock <- EuStockMarkets |>
#' #   ts_as_tibble() |>
#' #   group_by(stock_index = group)
#' # groups <- group_data(eu_stock) # Group information
#' # group_order <- unlist(groups$.rows) # Order of groups
#' # group_sizes <- lengths_(groups$.rows) # Group sizes
#' #
#' # # by-stock index lag
#' # lag2_(eu_stock$value, order = group_order, run_lengths = group_sizes)
#' #
#' # # Verifying this output is correct
#' # eu_stock |>
#' #   ungroup() |>
#' #   mutate(lag1 = lag_(value), .by = stock_index) |>
#' #   mutate(lag2 = lag2_(value, order = group_order, run_lengths = group_sizes)) |>
#' #   summarise(lags_are_equal = identical(lag1, lag2))
#'
#' # Let's compare this to data.table
#'
#' library(data.table)
#' setDTthreads(2)
#' dt <- data.table(x = 1:10^5,
#'                  g = sample.int(10^4, 10^5, TRUE))
#'
#' bench::mark(dt[, y := shift(x), by = g][][["y"]],
#'             grouped_lag(dt$x, g = dt$g),
#'             iterations = 10)
#' @rdname lag
#' @export
lag_ <- function(x, n = 1L, fill = NULL, set = FALSE, recursive = TRUE){
  .Call(`_cheapr_cpp_lag`, x, n, fill, set, recursive)
}
#' @rdname lag
#' @export
lag2_ <- function(x, n = 1L, order = NULL, run_lengths = NULL, fill = NULL, recursive = TRUE){
  .Call(`_cheapr_cpp_lag2`, x, n, order, run_lengths, fill, recursive)
}
# lag2_ <- function(x, n = 1L, order = if (recursive && is.data.frame(x)) seq_len(nrow(x)) else seq_along(x),
#                   run_lengths = if (recursive && is.data.frame(x)) nrow(x) else length(x),
#                   fill = NULL, recursive = TRUE){
#   .Call(`_cheapr_cpp_lag2`, x, n, order, run_lengths, fill, recursive)
# }

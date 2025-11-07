#' Extra utilities
#'
#' @param x A vector or data frame.
#' @param y A vector or data frame.
#' @param dups Should duplicates be kept? Default is `TRUE`.
#' @param value The column name to assign the values of a vector.
#' @param table See `?collapse::fmatch`
#' @param ... Further arguments passed onto `cut` or `set.seed`.
#' @param size See `?sample`.
#' @param replace See `?sample`.
#' @param prob See `?sample`.
#' @param n Number of scalar values (or `NA`) to insert
#' randomly into your vector.
#' @param prop Proportion of scalar values (or `NA`) values to insert
#' randomly into your vector.
#' @param na.rm Should `NA` values be ignored in `var_()` Default is
#' `TRUE`.
#' @param expr Expression that will be evaluated with a local seed that
#' is independent and has absolutely no effect on the global RNG state.
#' @param .seed A local seed to set which is only used inside
#' `with_local_seed()`. After the execution of the expression the original
#' seed is reset.
#' @param .envir Environment to evaluate expression.
#'
#' @returns
#' `intersect_()` returns a vector of common values between `x` and `y`. \cr
#' `setdiff_()` returns a vector of values in `x` but not `y`. \cr
#' `%in_%` and `%!in_%` both return a logical vector signifying if the values of
#' `x` exist or don't exist in `table` respectively. \cr
#' `sample_()` is an alternative to `sample()` that natively samples
#' data frame rows through `sset()`. It also does not have a special case when
#' `length(x)` is 1. \cr
#' `val_insert` inserts scalar values randomly into your vector.
#' Useful for replacing lots of data with a single value. \cr
#' `na_insert` inserts `NA` values randomly into your vector.
#' Useful for generating missing data. \cr
#' `var_` returns the variance of a numeric vector.
#' No coercion happens for integer vectors and so is very cheap. \cr
#' `rev_` is a much cheaper version of `rev()`. \cr
#' `with_local_seed` offers no speed improvements but is extremely handy
#' in executing random number based expressions like `rnorm()` without
#' affecting the global RNG state. It allows you to run these expressions in a
#' sort of independent 'container' and with an optional seed for that
#' 'container' for reproducibility.
#' The rationale for including this in 'cheapr' is that it can reduce the need
#' to set many seed values,
#' especially for multiple output comparisons of RNG expressions.
#' Another way of thinking about it is that `with_local_seed()` is a helper
#' that allows you to write reproducible code without side-effects, which
#' traditionally cannot be avoided when calling `set.seed()` directly.
#'
#' @examples
#' library(cheapr)
#'
#' # Using `with_local_seed()`
#'
#' # The below 2 statements are equivalent
#'
#' # Statement 1
#' set.seed(123456789)
#' res <- rnorm(10)
#'
#' # Statement 2
#' res2 <- with_local_seed(rnorm(10), .seed = 123456789)
#'
#' # They are the same
#' identical(res, res2)
#'
#' # As an example we can see that the RNG is unaffected by generating
#' # random uniform deviates in batches between calls to `with_local_seed()`
#' # and comparing to the first result
#'
#' set.seed(123456789)
#' batch1 <- rnorm(2)
#'
#' with_local_seed(runif(10))
#' batch2 <- rnorm(2)
#' with_local_seed(runif(10))
#' batch3 <- rnorm(1)
#' with_local_seed(runif(10))
#' batch4 <- rnorm(5)
#'
#' # Combining the batches produces the same result
#' # therefore `with_local_seed` did not interrupt the rng sequence
#' identical(c(batch1, batch2, batch3, batch4), res)
#'
#' # It can be useful in multiple comparisons
#' out1 <- with_local_seed(rnorm(5))
#' out2 <- with_local_seed(rnorm(5))
#' out3 <- with_local_seed(rnorm(5))
#'
#' identical(out1, out2)
#' identical(out1, out3)
#'
#' @rdname extras
#' @export
setdiff_ <- function(x, y, dups = TRUE){
  if (!dups){
    x <- unique_(x)
  }
  sset(x, which_not_in(x, y))
}
#' @rdname extras
#' @export
intersect_ <- function(x, y, dups = TRUE){
  if (!dups){
    x <- unique_(x)
  }
  sset(x, which_in(x, y))
}
#' @rdname extras
#' @export
`%in_%` <- function(x, table){
  collapse::fmatch(x, table, overid = 2L, nomatch = 0L) > 0L
}
#' @rdname extras
#' @export
`%!in_%` <- function(x, table){
  cpp_is_na(collapse::fmatch(x, table, overid = 2L, nomatch = NA_integer_))
}
cheapr_enframe <- function(x, name = "name", value = "value"){
  if (inherits(x, "data.frame")) {
    x <- unclass(x)
    attr(x, "row.names") <- NULL
  }
  x_nms <- names(x)
  x <- unname(x)
  if (is.null(x_nms)) {
    out <- list(x)
    names(out) <- value
  }
  else {
    out <- list(x_nms, x)
    names(out) <- c(name, value)
  }
  class(out) <- c("tbl_df", "tbl", "data.frame")
  attr(out, "row.names") <- .set_row_names(length(x))
  out
}
#' @rdname extras
#' @export
sample_ <- function(x, size = vector_length(x), replace = FALSE, prob = NULL){
  sset(x, sample.int(vector_length(x), size, replace, prob))
}
#' @rdname extras
#' @export
val_insert <- function(x, value, n = NULL, prop = NULL){
  if (!is.null(n) && !is.null(prop)) {
    stop("either n or prop must be supplied")
  }
  if (!is.null(n)){
    x[sample.int(length(x), size = n, replace = FALSE)] <- value
  }
  if (!is.null(prop)) {
    x[sample.int(length(x), size = floor(prop * length(x)),
                 replace = FALSE)] <- value
  }
  x
}
#' @rdname extras
#' @export
na_insert <- function(x, n = NULL, prop = NULL){
  val_insert(x, value = NA, n = n, prop = prop)
}
#' @rdname extras
#' @export
vector_length <- cpp_vector_length
#' @rdname extras
#' @export
cheapr_var <- function(x, na.rm = TRUE){
  y <- as_numeric(x)
  if (na.rm){
    N <- (length(y) - na_count(y)) - 1
  } else {
    N <- length(y) - 1
  }
  if (N < 1){
    return(NA_real_)
  } else {
    mu <- collapse::fmean(y, na.rm = na.rm)
    # mu <- sum(y, na.rm = na.rm) / (N + 1)
    if (is.na(mu)){
      NA_real_
    } else {
      var_sum_squared_diff(y, as.double(mu)) /  N
    }
  }
}
#' @rdname extras
#' @export
cheapr_rev <- function(x){
  .Call(`_cheapr_cpp_sset`, x, .Call(`_cheapr_cpp_vector_length`, x):0, TRUE)
}
#' @rdname extras
#' @export
cheapr_sd <- function(x, na.rm = TRUE){
  sqrt(cheapr_var(x, na.rm = na.rm))
}

#' @rdname extras
#' @export
rev_ <- cheapr_rev
#' @rdname extras
#' @export
sd_ <- cheapr_sd
#' @rdname extras
#' @export
var_ <- cheapr_var
#' @rdname extras
#' @export
with_local_seed <- function (expr, .seed = NULL, .envir = environment(), ...){
  old <- globalenv()[[".Random.seed"]]
  if (is.null(old)){
    set.seed(NULL)
    old <- globalenv()[[".Random.seed"]]
  }
  if (!is.null(.seed)){
    set.seed(.seed, ...)
  }
  on.exit({
    assign(".Random.seed", old, envir = globalenv())
  }, add = TRUE)
  eval(expr, envir = .envir)
}

# is_duplicate <- function(x, .all = FALSE){
#   groups <- collapse::group(x, starts = !.all, group.sizes = TRUE)
#   sizes <- attr(groups, "group.sizes", TRUE)
#   starts <- attr(groups, "starts", TRUE)
#   out <- (sizes > 1L)[groups]
#   out[starts] <- FALSE
#   out
# }
# duplicates <- function(x, .all = FALSE, .count = FALSE){
#   groups <- collapse::group(x, starts = !.all, group.sizes = TRUE)
#   sizes <- attr(groups, "group.sizes")
#   starts <- attr(groups, "starts")
#   dup <- (sizes > 1L)[groups]
#   dup[starts] <- FALSE
#   which_dup <- which_(dup)
#   out <- sset(x, which_dup)
#
#   # Adjust group sizes as they reflect the dup count + 1
#
#   if (.count){
#     sizes <- sizes[groups]
#     if (!.all && NROW(out) > 0){
#       set_subtract(sizes, 1L)
#       which_zero <- which_val(sizes, 0L)
#       collapse::setv(
#         sizes,
#         which_zero,
#         1L,
#         vind1 = TRUE
#       )
#     }
#     cpp_set_add_attr(out, "n_dupes", sset(sizes, which_dup))
#   }
#   out
# }

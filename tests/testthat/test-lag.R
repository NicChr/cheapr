base_lag <- function(x, lag = 1L, check = TRUE){
  lagged_indices <- seq_along(x) - lag
  lagged_indices[lagged_indices < 1L] <- NA_integer_
  x[lagged_indices]
}

test_that("lags and leads", {
  base_lag <- function(x, lag = 1L, check = TRUE){
    lagged_indices <- seq_along(x) - lag
    lagged_indices[lagged_indices < 1L] <- NA_integer_
    x[lagged_indices]
  }
  set.seed(876123)
  a <- rnorm(10^5)
  b <- sample(-12:123, 10^5, TRUE)
  c <- sample(do.call(paste0, expand.grid(letters, letters, letters)), 10^5, TRUE)
  d <- complex(na_insert(rnorm(10^5, 10^3)),
               na_insert(rnorm(10^5, 10^3)))
  e <- vapply(sample(letters, 10^5, TRUE), charToRaw, raw(1))
  a <- na_insert(a, 10^3)
  names(a) <- sample(letters, length(a), TRUE)
  b <- na_insert(b, 10^3)
  names(b) <- sample(letters, length(b), TRUE)
  c <- na_insert(c, 10^3)
  names(c) <- sample(letters, length(c), TRUE)
  names(d) <- sample(letters, length(d), TRUE)
  f <- as.list(b)
  g <- a > 0

  expect_identical(
   lag_(a, 0), a
  )
  expect_identical(
    lag_(a, 1), base_lag(a, 1)
  )
  expect_identical(
    lag_(a, -1), base_lag(a, -1)
  )
  expect_identical(
    lag_(a, 3), base_lag(a, 3)
  )
  expect_identical(
    lag_(a, -3), base_lag(a, -3)
  )
  expect_identical(
    lag_(a, 10^6), base_lag(a, 10^6)
  )
  expect_identical(
    lag_(a, -10^6), base_lag(a, -10^6)
  )

  expect_identical(
    lag_(b, 0), b
  )
  expect_identical(
    lag_(b, 1), base_lag(b, 1)
  )
  expect_identical(
    lag_(b, -1), base_lag(b, -1)
  )
  expect_identical(
    lag_(b, 3), base_lag(b, 3)
  )
  expect_identical(
    lag_(b, -3), base_lag(b, -3)
  )
  expect_identical(
    lag_(b, 10^6), base_lag(b, 10^6)
  )
  expect_identical(
    lag_(b, -10^6), base_lag(b, -10^6)
  )


  expect_identical(
    lag_(c, 0), c
  )
  expect_identical(
    lag_(c, 1), base_lag(c, 1)
  )
  expect_identical(
    lag_(c, -1), base_lag(c, -1)
  )
  expect_identical(
    lag_(c, 3), base_lag(c, 3)
  )
  expect_identical(
    lag_(c, -3), base_lag(c, -3)
  )
  expect_identical(
    lag_(c, 10^6), base_lag(c, 10^6)
  )
  expect_identical(
    lag_(c, -10^6), base_lag(c, -10^6)
  )


  expect_identical(
    lag_(d, 0), d
  )
  expect_identical(
    lag_(d, 1), base_lag(d, 1)
  )
  expect_identical(
    lag_(d, -1), base_lag(d, -1)
  )
  expect_identical(
    lag_(d, 3), base_lag(d, 3)
  )
  expect_identical(
    lag_(d, -3), base_lag(d, -3)
  )
  expect_identical(
    lag_(d, 10^6), base_lag(d, 10^6)
  )
  expect_identical(
    lag_(d, -10^6), base_lag(d, -10^6)
  )


  expect_identical(
    lag_(e, 0), e
  )
  expect_identical(
    lag_(e, 1), base_lag(e, 1)
  )
  expect_identical(
    lag_(e, -1), base_lag(e, -1)
  )
  expect_identical(
    lag_(e, 3), base_lag(e, 3)
  )
  expect_identical(
    lag_(e, -3), base_lag(e, -3)
  )
  expect_identical(
    lag_(e, 10^6), base_lag(e, 10^6)
  )
  expect_identical(
    lag_(e, -10^6), base_lag(e, -10^6)
  )

  expect_identical(
    lag_(f, 0), f
  )
  expect_identical(
    lag_(f, 1, recursive = FALSE), base_lag(f, 1)
  )
  expect_identical(
    lag_(f, -1, recursive = FALSE), base_lag(f, -1)
  )
  expect_identical(
    lag_(f, 3, recursive = FALSE), base_lag(f, 3)
  )
  expect_identical(
    lag_(f, -3, recursive = FALSE), base_lag(f, -3)
  )
  expect_identical(
    lag_(f, 10^6, recursive = FALSE), base_lag(f, 10^6)
  )
  expect_identical(
    lag_(f, -10^6, recursive = FALSE), base_lag(f, -10^6)
  )

  expect_identical(
    lag_(g, 0), g
  )
  expect_identical(
    lag_(g, 1, recursive = FALSE), base_lag(g, 1)
  )
  expect_identical(
    lag_(g, -1, recursive = FALSE), base_lag(g, -1)
  )
  expect_identical(
    lag_(g, 3, recursive = FALSE), base_lag(g, 3)
  )
  expect_identical(
    lag_(g, -3, recursive = FALSE), base_lag(g, -3)
  )
  expect_identical(
    lag_(g, 10^6, recursive = FALSE), base_lag(g, 10^6)
  )
  expect_identical(
    lag_(g, -10^6, recursive = FALSE), base_lag(g, -10^6)
  )

  expect_identical(
    lag_(iris, 7),
    as.data.frame(lapply(iris, base_lag, 7)),
  )
})

test_that("lags and leads with lag2_", {
  base_lag <- function(x, lag = 1L, check = TRUE){
    lagged_indices <- seq_along(x) - lag
    lagged_indices[lagged_indices < 1L] <- NA_integer_
    x[lagged_indices]
  }
  set.seed(876123)
  a <- rnorm(10^5)
  b <- sample(-12:123, 10^5, TRUE)
  c <- sample(do.call(paste0, expand.grid(letters, letters, letters)), 10^5, TRUE)
  d <- complex(na_insert(rnorm(10^5, 10^3)),
               na_insert(rnorm(10^5, 10^3)))
  e <- vapply(sample(letters, 10^5, TRUE), charToRaw, raw(1))
  a <- na_insert(a, 10^3)
  names(a) <- sample(letters, length(a), TRUE)
  b <- na_insert(b, 10^3)
  names(b) <- sample(letters, length(b), TRUE)
  c <- na_insert(c, 10^3)
  names(c) <- sample(letters, length(c), TRUE)
  names(d) <- sample(letters, length(d), TRUE)
  f <- as.list(b)
  g <- a > 0

  expect_identical(
    lag2_(a, 0), a
  )
  expect_identical(
    lag2_(a, 1), base_lag(a, 1)
  )
  expect_identical(
    lag2_(a, -1), base_lag(a, -1)
  )
  expect_identical(
    lag2_(a, 3), base_lag(a, 3)
  )
  expect_identical(
    lag2_(a, -3), base_lag(a, -3)
  )
  expect_identical(
    lag2_(a, 10^6), base_lag(a, 10^6)
  )
  expect_identical(
    lag2_(a, -10^6), base_lag(a, -10^6)
  )

  expect_identical(
    lag2_(b, 0), b
  )
  expect_identical(
    lag2_(b, 1), base_lag(b, 1)
  )
  expect_identical(
    lag2_(b, -1), base_lag(b, -1)
  )
  expect_identical(
    lag2_(b, 3), base_lag(b, 3)
  )
  expect_identical(
    lag2_(b, -3), base_lag(b, -3)
  )
  expect_identical(
    lag2_(b, 10^6), base_lag(b, 10^6)
  )
  expect_identical(
    lag2_(b, -10^6), base_lag(b, -10^6)
  )


  expect_identical(
    lag2_(c, 0), c
  )
  expect_identical(
    lag2_(c, 1), base_lag(c, 1)
  )
  expect_identical(
    lag2_(c, -1), base_lag(c, -1)
  )
  expect_identical(
    lag2_(c, 3), base_lag(c, 3)
  )
  expect_identical(
    lag2_(c, -3), base_lag(c, -3)
  )
  expect_identical(
    lag2_(c, 10^6), base_lag(c, 10^6)
  )
  expect_identical(
    lag2_(c, -10^6), base_lag(c, -10^6)
  )


  expect_identical(
    lag2_(d, 0), d
  )
  expect_identical(
    lag2_(d, 1), base_lag(d, 1)
  )
  expect_identical(
    lag2_(d, -1), base_lag(d, -1)
  )
  expect_identical(
    lag2_(d, 3), base_lag(d, 3)
  )
  expect_identical(
    lag2_(d, -3), base_lag(d, -3)
  )
  expect_identical(
    lag2_(d, 10^6), base_lag(d, 10^6)
  )
  expect_identical(
    lag2_(d, -10^6), base_lag(d, -10^6)
  )


  expect_identical(
    lag2_(e, 0), e
  )
  expect_identical(
    lag2_(e, 1), base_lag(e, 1)
  )
  expect_identical(
    lag2_(e, -1), base_lag(e, -1)
  )
  expect_identical(
    lag2_(e, 3), base_lag(e, 3)
  )
  expect_identical(
    lag2_(e, -3), base_lag(e, -3)
  )
  expect_identical(
    lag2_(e, 10^6), base_lag(e, 10^6)
  )
  expect_identical(
    lag2_(e, -10^6), base_lag(e, -10^6)
  )

  expect_identical(
    lag2_(f, 0), f
  )
  expect_identical(
    lag2_(f, 1, recursive = FALSE), base_lag(f, 1)
  )
  expect_identical(
    lag2_(f, -1, recursive = FALSE), base_lag(f, -1)
  )
  expect_identical(
    lag2_(f, 3, recursive = FALSE), base_lag(f, 3)
  )
  expect_identical(
    lag2_(f, -3, recursive = FALSE), base_lag(f, -3)
  )
  expect_identical(
    lag2_(f, 10^6, recursive = FALSE), base_lag(f, 10^6)
  )
  expect_identical(
    lag2_(f, -10^6, recursive = FALSE), base_lag(f, -10^6)
  )

  expect_identical(
    lag2_(g, 0), g
  )
  expect_identical(
    lag2_(g, 1, recursive = FALSE), base_lag(g, 1)
  )
  expect_identical(
    lag2_(g, -1, recursive = FALSE), base_lag(g, -1)
  )
  expect_identical(
    lag2_(g, 3, recursive = FALSE), base_lag(g, 3)
  )
  expect_identical(
    lag2_(g, -3, recursive = FALSE), base_lag(g, -3)
  )
  expect_identical(
    lag2_(g, 10^6, recursive = FALSE), base_lag(g, 10^6)
  )
  expect_identical(
    lag2_(g, -10^6, recursive = FALSE), base_lag(g, -10^6)
  )

  expect_identical(
    lag2_(iris, 7),
    as.data.frame(lapply(iris, base_lag, 7)),
  )
})

test_that("lags and lead with set = TRUE", {
  set.seed(876123)
  a <- rnorm(10^5)
  b <- sample(-12:123, 10^5, TRUE)
  c <- sample(do.call(paste0, expand.grid(letters, letters, letters)), 10^5, TRUE)
  d <- complex(na_insert(rnorm(10^5, 10^3)),
               na_insert(rnorm(10^5, 10^3)))
  e <- vapply(sample(letters, 10^5, TRUE), charToRaw, raw(1))
  a <- na_insert(a, 10^3)
  names(a) <- sample(letters, length(a), TRUE)
  b <- na_insert(b, 10^3)
  names(b) <- sample(letters, length(b), TRUE)
  c <- na_insert(c, 10^3)
  names(c) <- sample(letters, length(c), TRUE)
  names(d) <- sample(letters, length(d), TRUE)
  f <- as.list(b)
  g <- a > 0

  set_lag <- function(x, ...){
    lag_(x, ..., set = TRUE)
  }

  expect_identical(
    set_lag(deep_copy(a), 0), a
  )
  expect_identical(
    set_lag(deep_copy(a), 1), base_lag(a, 1)
  )
  expect_identical(
    set_lag(deep_copy(a), -1), base_lag(a, -1)
  )
  expect_identical(
    set_lag(deep_copy(a), 3), base_lag(a, 3)
  )
  expect_identical(
    set_lag(deep_copy(a), -3), base_lag(a, -3)
  )
  expect_identical(
    set_lag(deep_copy(a), 10^6), base_lag(a, 10^6)
  )
  expect_identical(
    set_lag(deep_copy(a), -10^6), base_lag(a, -10^6)
  )

  expect_identical(
    set_lag(deep_copy(b), 0), b
  )
  expect_identical(
    set_lag(deep_copy(b), 1), base_lag(b, 1)
  )
  expect_identical(
    set_lag(deep_copy(b), -1), base_lag(b, -1)
  )
  expect_identical(
    set_lag(deep_copy(b), 3), base_lag(b, 3)
  )
  expect_identical(
    set_lag(deep_copy(b), -3), base_lag(b, -3)
  )
  expect_identical(
    set_lag(deep_copy(b), 10^6), base_lag(b, 10^6)
  )
  expect_identical(
    set_lag(deep_copy(b), -10^6), base_lag(b, -10^6)
  )


  expect_identical(
    set_lag(deep_copy(c), 0), c
  )
  expect_identical(
    set_lag(deep_copy(c), 1), base_lag(c, 1)
  )
  expect_identical(
    set_lag(deep_copy(c), -1), base_lag(c, -1)
  )
  expect_identical(
    set_lag(deep_copy(c), 3), base_lag(c, 3)
  )
  expect_identical(
    set_lag(deep_copy(c), -3), base_lag(c, -3)
  )
  expect_identical(
    set_lag(deep_copy(c), 10^6), base_lag(c, 10^6)
  )
  expect_identical(
    set_lag(deep_copy(c), -10^6), base_lag(c, -10^6)
  )


  expect_identical(
    set_lag(deep_copy(d), 0), d
  )
  expect_identical(
    set_lag(deep_copy(d), 1), base_lag(d, 1)
  )
  expect_identical(
    set_lag(deep_copy(d), -1), base_lag(d, -1)
  )
  expect_identical(
    set_lag(deep_copy(d), 3), base_lag(d, 3)
  )
  expect_identical(
    set_lag(deep_copy(d), -3), base_lag(d, -3)
  )
  expect_identical(
    set_lag(deep_copy(d), 10^6), base_lag(d, 10^6)
  )
  expect_identical(
    set_lag(deep_copy(d), -10^6), base_lag(d, -10^6)
  )


  expect_identical(
    set_lag(deep_copy(e), 0), e
  )
  expect_identical(
    set_lag(deep_copy(e), 1), base_lag(e, 1)
  )
  expect_identical(
    set_lag(deep_copy(e), -1), base_lag(e, -1)
  )
  expect_identical(
    set_lag(deep_copy(e), 3), base_lag(e, 3)
  )
  expect_identical(
    set_lag(deep_copy(e), -3), base_lag(e, -3)
  )
  expect_identical(
    set_lag(deep_copy(e), 10^6), base_lag(e, 10^6)
  )
  expect_identical(
    set_lag(deep_copy(e), -10^6), base_lag(e, -10^6)
  )

  expect_identical(
    set_lag(deep_copy(f), 0), f
  )
  expect_identical(
    set_lag(deep_copy(f), 1, recursive = FALSE), base_lag(f, 1)
  )
  expect_identical(
    set_lag(deep_copy(f), -1, recursive = FALSE), base_lag(f, -1)
  )
  expect_identical(
    set_lag(deep_copy(f), 3, recursive = FALSE), base_lag(f, 3)
  )
  expect_identical(
    set_lag(deep_copy(f), -3, recursive = FALSE), base_lag(f, -3)
  )
  expect_identical(
    set_lag(deep_copy(f), 10^6, recursive = FALSE), base_lag(f, 10^6)
  )
  expect_identical(
    set_lag(deep_copy(f), -10^6, recursive = FALSE), base_lag(f, -10^6)
  )

  expect_identical(
    set_lag(deep_copy(g), 0), g
  )
  expect_identical(
    set_lag(deep_copy(g), 1, recursive = FALSE), base_lag(g, 1)
  )
  expect_identical(
    set_lag(deep_copy(g), -1, recursive = FALSE), base_lag(g, -1)
  )
  expect_identical(
    set_lag(deep_copy(g), 3, recursive = FALSE), base_lag(g, 3)
  )
  expect_identical(
    set_lag(deep_copy(g), -3, recursive = FALSE), base_lag(g, -3)
  )
  expect_identical(
    set_lag(deep_copy(g), 10^6, recursive = FALSE), base_lag(g, 10^6)
  )
  expect_identical(
    set_lag(deep_copy(g), -10^6, recursive = FALSE), base_lag(g, -10^6)
  )

  expect_identical(
    set_lag(deep_copy(iris), 7),
    as.data.frame(lapply(iris, base_lag, 7)),
  )
})

test_that("Dynamic lags by-group", {
  set.seed(1239)
  df <- data.frame(x = sample.int(5, 20, TRUE),
                   g = sample.int(3, 20, TRUE),
                   lags = sample(c(0, 1, 2), 20, TRUE))

  o <- order(df$g)
  rls <- as.integer(table(df$g))

  # Somewhat ugly by-group calculation
  # order(order(x)) will return sort(x) back to its original order
  res <- unname(
    do.call(
      c,
      lapply(split(df, df$g),
             function(x) base_lag(x$x, x$lags))
    )
  )[order(o)]

res2 <- lag2_(df$x, order = o, run_lengths = rls, n = df$lags)

expect_identical(res, res2)
})

test_that("oob lag", {
  expect_identical(lag_(1:10, 100), rep(NA_integer_, 10))
  expect_identical(lag_(1:10, -100), rep(NA_integer_, 10))
  expect_identical(lag_(1:10, 100, fill = 99), rep(99L, 10))
  expect_identical(lag_(1:10, -100, fill = 99), rep(99L, 10))
})

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
  d <- complex(fill_with_na(rnorm(10^5, 10^3)),
               fill_with_na(rnorm(10^5, 10^3)))
  e <- vapply(sample(letters, 10^5, TRUE), charToRaw, raw(1))
  a <- fill_with_na(a, 10^3)
  names(a) <- sample(letters, length(a), TRUE)
  b <- fill_with_na(b, 10^3)
  names(b) <- sample(letters, length(b), TRUE)
  c <- fill_with_na(c, 10^3)
  names(c) <- sample(letters, length(c), TRUE)
  names(d) <- sample(letters, length(d), TRUE)
  f <- as.list(b)

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
    lag_(iris, 7),
    as.data.frame(lapply(iris, base_lag, 7)),
  )
})

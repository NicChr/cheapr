test_that("GCD AND SCM scalars", {
  expect_identical(gcd(1), 1)
  expect_identical(gcd(1L), 1L)
  expect_identical(gcd(0), 0)
  expect_identical(gcd(0L), 0L)
  expect_identical(gcd(NaN), NaN)
  expect_identical(gcd(Inf), Inf)
  expect_identical(gcd(-Inf), Inf)
  expect_identical(gcd(NA_integer_, na_rm = TRUE), NA_integer_)
  expect_identical(gcd(NA_integer_, na_rm = FALSE), NA_integer_)
  expect_identical(gcd(NA_real_, na_rm = TRUE), NA_real_)
  expect_identical(gcd(NA_real_, na_rm = FALSE), NA_real_)

  expect_identical(scm(1), 1)
  expect_identical(scm(1L), 1L)
  expect_identical(scm(0), 0)
  expect_identical(scm(0L), 0L)
  expect_identical(scm(NaN), NaN)
  expect_identical(scm(Inf), Inf)
  expect_identical(scm(-Inf), -Inf)
  expect_identical(scm(NA_integer_, na_rm = TRUE), NA_integer_)
  expect_identical(scm(NA_integer_, na_rm = FALSE), NA_integer_)
  expect_identical(scm(NA_real_, na_rm = TRUE), NA_real_)
  expect_identical(scm(NA_real_, na_rm = FALSE), NA_real_)
})

test_that("zero-length vectors", {
  expect_identical(gcd(integer()), integer())
  expect_identical(gcd(integer(), na_rm = FALSE), integer())
  expect_identical(gcd(numeric()), numeric())
  expect_identical(gcd(numeric(), na_rm = FALSE), numeric())

  expect_identical(scm(integer()), integer())
  expect_identical(scm(integer(), na_rm = FALSE), integer())
  expect_identical(scm(numeric()), numeric())
  expect_identical(scm(numeric(), na_rm = FALSE), numeric())
})

test_that("Small vectors (integers)", {
  x <- c(NA_integer_, 9L, 0L, 0L, 3L, 21L)
  y <- c(9L, NA_integer_, NA_integer_, 0L, 0L, 0L, 3L, 21L)
  expect_identical(gcd(x, na_rm = TRUE), 3L)
  expect_identical(gcd(x, na_rm = FALSE), NA_integer_)
  expect_identical(gcd(x, tol = 0, na_rm = TRUE), 3L)
  expect_identical(gcd(x, tol = 0, na_rm = FALSE), NA_integer_)
  expect_identical(gcd(x, tol = 0, na_rm = TRUE, round = TRUE), 3L)

  expect_identical(gcd(y, na_rm = TRUE), 3L)
  expect_identical(gcd(y, na_rm = FALSE), NA_integer_)
  expect_identical(gcd(y, tol = 0, na_rm = TRUE), 3L)
  expect_identical(gcd(y, tol = 0, na_rm = FALSE), NA_integer_)
  expect_identical(gcd(y, tol = 0, na_rm = TRUE, round = TRUE), 3L)

  expect_identical(gcd(as.integer(c(NA, 14, 21)), na_rm = TRUE), 7L)
  expect_identical(gcd(as.integer(c(21, NA, 14)), na_rm = TRUE), 7L)
  expect_identical(gcd(as.integer(c(14, 21, NA)), na_rm = TRUE), 7L)
  expect_identical(gcd(as.integer(c(NA, 21, 14)), na_rm = FALSE), NA_integer_)
  expect_identical(gcd(as.integer(c(21, NA, 14)), na_rm = FALSE), NA_integer_)
  expect_identical(gcd(as.integer(c(14, 21, NA)), na_rm = FALSE), NA_integer_)

  expect_identical(gcd(as.integer(c(NA, NA, 14, NA, 21, NA, NA)), na_rm = TRUE), 7L)
  expect_identical(gcd(as.integer(c(NA, NA, 14, NA, 21, NA, NA)), na_rm = FALSE), NA_integer_)

  expect_identical(scm(as.integer(c(NA, 2, 3)), na_rm = TRUE), 6L)
  expect_identical(scm(as.integer(c(2, NA, 3)), na_rm = TRUE), 6L)
  expect_identical(scm(as.integer(c(2, 3, NA)), na_rm = TRUE), 6L)
  expect_identical(scm(as.integer(c(NA, 2, 3)), na_rm = FALSE), NA_integer_)
  expect_identical(scm(as.integer(c(2, NA, 3)), na_rm = FALSE), NA_integer_)
  expect_identical(scm(as.integer(c(2, 3, NA)), na_rm = FALSE), NA_integer_)

  expect_identical(scm(as.integer(c(NA, NA, 2, NA, 3, NA, NA)), na_rm = TRUE), 6L)
  expect_identical(scm(as.integer(c(NA, NA, 2, NA, 3, NA, NA)), na_rm = FALSE), NA_integer_)

  expect_identical(gcd(c(0, Inf)), Inf)
  expect_identical(gcd(c(Inf, 0)), Inf)
  expect_identical(gcd(c(0, -Inf)), Inf)
  expect_identical(gcd(c(-Inf, 0)), Inf)
  expect_identical(scm(c(0, Inf)), NaN)
  expect_identical(scm(c(Inf, 0)), NaN)
  expect_identical(scm(c(0, -Inf)), NaN)
  expect_identical(scm(c(-Inf, 0)), NaN)
})

test_that("Small vectors (doubles)", {
  x <- c(NA_real_, 9, 0, 0, 3, 21)
  y <- c(9, NA_real_, NA_real_, 0, 0, 0, NaN, NaN, NaN, 3, Inf, -Inf, 21)
  expect_identical(gcd(x, na_rm = TRUE), 3)
  expect_identical(gcd(x, na_rm = FALSE), NA_real_)
  expect_identical(gcd(x, tol = 0, na_rm = TRUE), 3)
  expect_identical(gcd(x, tol = 0, na_rm = FALSE), NA_real_)
  expect_identical(gcd(x, tol = 0, na_rm = TRUE, round = TRUE), 3)

  expect_identical(gcd(y, na_rm = TRUE), 3)
  expect_identical(gcd(y, na_rm = FALSE), NA_real_)
  expect_identical(gcd(y, tol = 0, na_rm = TRUE), 3)
  expect_identical(gcd(y, tol = 0, na_rm = FALSE), NA_real_)
  expect_identical(gcd(y, tol = 0, na_rm = TRUE, round = TRUE), 3)

  expect_identical(gcd(as.double(c(NA, 14, 21)), na_rm = TRUE), 7)
  expect_identical(gcd(as.double(c(21, NA, 14)), na_rm = TRUE), 7)
  expect_identical(gcd(as.double(c(14, 21, NA)), na_rm = TRUE), 7)
  expect_identical(gcd(as.double(c(NA, 21, 14)), na_rm = FALSE), NA_real_)
  expect_identical(gcd(as.double(c(21, NA, 14)), na_rm = FALSE), NA_real_)
  expect_identical(gcd(as.double(c(14, 21, NA)), na_rm = FALSE), NA_real_)

  expect_identical(gcd(as.double(c(NA, NA, 14, NA, 21, NA, NA)), na_rm = TRUE), 7)
  expect_identical(gcd(as.double(c(NA, NA, 14, NA, 21, NA, NA)), na_rm = FALSE), NA_real_)

  expect_identical(scm(as.double(c(NA, 2, 3)), na_rm = TRUE), 6)
  expect_identical(scm(as.double(c(2, NA, 3)), na_rm = TRUE), 6)
  expect_identical(scm(as.double(c(2, 3, NA)), na_rm = TRUE), 6)
  expect_identical(scm(as.double(c(NA, 2, 3)), na_rm = FALSE), NA_real_)
  expect_identical(scm(as.double(c(2, NA, 3)), na_rm = FALSE), NA_real_)
  expect_identical(scm(as.double(c(2, 3, NA)), na_rm = FALSE), NA_real_)

  expect_identical(scm(as.double(c(NA, NA, 2, NA, 3, NA, NA)), na_rm = TRUE), 6)
  expect_identical(scm(as.double(c(NA, NA, 2, NA, 3, NA, NA)), na_rm = FALSE), NA_real_)
})

test_that("Small vectors (doubles)", {
  set.seed(94773)
  x <- rnorm(100)
  x <- abs(x)
  tol <- sqrt(.Machine$double.eps)
  expect_equal(gcd(x, tol = tol, break_early = TRUE), tol)
  expect_equal(gcd(x, tol = tol, break_early = FALSE), tol)
  expect_equal(gcd(x, tol = tol, round = TRUE), round(tol, 9))
  expect_equal(gcd(x, tol = tol, round = TRUE, break_early = TRUE),
               round(tol, 9))
  expect_true(gcd(x, tol = 0) < tol/10)

  x <- as.double(time(AirPassengers))

  expect_equal(gcd(x, tol = tol), 1/12)
  expect_false(isTRUE(all.equal(gcd(x, tol = 0), 1/12)))
  expect_false(isTRUE(all.equal(gcd(diff(x), tol = 0), 1/12)))
  expect_equal(gcd(x, round = TRUE, tol = tol), 1/12)
  expect_equal(gcd(diff(x), tol = tol), 1/12)
  expect_equal(gcd(diff(x), round = TRUE, tol = tol), 1/12)
})

test_that("signs", {
  a <- c(0L, 0L, seq.int(25L, 0L, by = -5L))
  b <- c(0, 0, seq.int(25, 0, by = -5))
  c <- -a
  d <- -b
  expect_identical(gcd(a), 5L)
  expect_identical(gcd(b), 5)
  expect_identical(gcd(c), 5L)
  expect_identical(gcd(d), 5)
})

test_that("Overflow", {
  expect_equal(scm(1:30), 2329089562800)
  expect_equal(scm(-(1:30)), 2329089562800)
  expect_equal(scm(1:50), 1.85233763824414e+28)
  expect_equal(scm(-(1:50)), 1.85233763824414e+28)
  expect_equal(scm(as.double(1:50)), 18523376382441352270484866200)
  expect_equal(scm(-as.double(1:50)), 18523376382441352270484866200)
})

test_that("Binary gcd and scm", {
  expect_identical(gcd2(10, NA_real_), 10)
  expect_identical(gcd2(10, NA, na_rm = FALSE), NA_real_)
  expect_identical(gcd2(10L, NA), 10L)
  expect_identical(gcd2(10L, NA_integer_, na_rm = FALSE), NA_integer_)
  expect_identical(gcd2(0, 0), 0)
  expect_identical(gcd2(0L, 0L), 0L)
  expect_identical(gcd2(3, 0), 3)
  expect_identical(gcd2(0, 3), 3)
  expect_identical(gcd2(3L, 0L), 3L)
  expect_identical(gcd2(0L, 3L), 3L)
  expect_identical(gcd2(1, 0), 1)
  expect_identical(gcd2(0, 1), 1)
  expect_identical(gcd2(1L, 0L), 1L)
  expect_identical(gcd2(0L, 1L), 1L)

  expect_identical(gcd2(seq(0, 25, 5), 5), rep(5, 6))
  expect_identical(gcd2(seq(0L, 25L, 5L), 5L), rep(5L, 6))

  expect_identical(scm2(0, 0), 0)
  expect_identical(scm2(0L, 0L), 0L)
  expect_identical(scm2(3, 0), 0)
  expect_identical(scm2(0, 3), 0)
  expect_identical(scm2(3L, 0L), 0L)
  expect_identical(scm2(0L, 3L), 0L)

  expect_identical(scm2(seq(0, 25, 0.5), 0.5), seq(0, 25, 0.5))
  expect_identical(scm2(seq(0L, 25L, 5L), 5L), seq(0L, 25L, 5L))

  expect_identical(scm2(3L, 1234567891L), NA_integer_)
  expect_identical(scm2(3, 1234567891), 3703703673)

  expect_identical(scm2(Inf, 0), NaN)
  expect_identical(scm2(Inf, 2), Inf)
  expect_identical(scm2(Inf, 2), Inf)
  expect_identical(gcd2(Inf, 2), 2)
  expect_identical(gcd2(Inf, 0), Inf) # Not sure about this
})

test_that("breaks", {
  x <- 1:20

  b1 <- function(x, n){
   get_breaks(x, n, expand_min = FALSE, expand_max = FALSE)
  }

  b2 <- function(x, n){
    get_breaks(x, n, expand_min = FALSE, expand_max = TRUE)
  }

  b3 <- function(x, n){
    get_breaks(x, n, expand_min = FALSE, expand_max = FALSE, pretty = FALSE)
  }

  expect_equal(b1(x, 5), seq(0, 20, 5))
  expect_equal(b2(x, 5), seq(0, 25, 5))

  expect_equal(b1(1:10, 5), seq(1, 9, 2))
  expect_equal(b2(1:10, 5), seq(1, 11, 2))

  expect_equal(b1(1:10, 20), seq(1, 10, 0.5))
  expect_equal(b2(1:10, 20), seq(1, 10.5, 0.5))

  expect_equal(b3(1:10, 20), seq(1, 10, 9/20))
  expect_equal(b3(1:10, 10), seq(1, 10, 9/10))
  expect_equal(b3(1:17, 13), seq(1, 17, 16/13))

  expect_equal(
    get_breaks(0:10, pretty = FALSE, expand_min = TRUE, expand_max = TRUE,
               n = 7),
    seq(0 - (10/7), 10 + (10/7), 10/7)
  )

  expect_equal(
    b1(with_local_seed(rnorm(100), .seed = 42), 13),
    seq(-3, 2, 0.5)
  )


  expect_equal(
    b1(seq(-0.034, 123, 0.5), 10),
    seq(-30, 120, 30)
  )

  expect_equal(get_breaks(1:30, 3), seq(0, 40, 10))
  expect_equal(get_breaks(1:30, 5), seq(0, 40, 10))

  # Zero range

  expect_equal(b1(123, 7), seq(123 - 123/1000, 123 + 123/1000, length.out = 8))

  # Lots of breaks

  expect_equal(
    b2(with_local_seed(rnorm(100), .seed = 42), 123456),
    seq(-2.99310, 2.28665, 0.00005)
  )

  # Floating point precision example

  expect_equal(
    get_breaks(with_local_seed(rnorm(100, sd = 0.01), .seed = 3), 10),
    seq(-0.025, 0.02, 0.005)
  )

})

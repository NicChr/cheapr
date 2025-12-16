options(cheapr.cores = 2)

test_that("binning", {
  set.seed(42)
  x <- sample(-10:11, 100, TRUE)
  breaks <- seq(0L, 7L, by = 1L)

  .bin <- function(x, breaks, ...){
    breaks[.bincode(x, breaks, ...)]
  }

  expect_equal(
    .bincode(x, breaks),
    bin(x, breaks, left_closed = FALSE)
  )
  expect_equal(
    .bincode(x, breaks, right = FALSE),
    bin(x, breaks, left_closed = TRUE)
  )
  expect_equal(
    .bincode(x, breaks, include.lowest = TRUE),
    bin(x, breaks, include_endpoint = TRUE, left_closed = FALSE)
  )
  expect_equal(
    .bincode(x, breaks, right = TRUE, include.lowest = TRUE),
    bin(x, breaks, left_closed = FALSE, include_endpoint = TRUE)
  )
  expect_equal(
    .bincode(x, breaks, right = TRUE, include.lowest = FALSE),
    bin(x, breaks, left_closed = FALSE, include_endpoint = FALSE)
  )

  breaks <- seq(0, max(x), by = 0.5)

  expect_equal(
    .bincode(x, breaks),
    bin(x, breaks, left_closed = FALSE)
  )
  expect_equal(
    .bincode(x, breaks, right = FALSE),
    bin(x, breaks, left_closed = TRUE)
  )
  expect_equal(
    .bincode(x, breaks, include.lowest = TRUE),
    bin(x, breaks, include_endpoint = TRUE, left_closed = FALSE)
  )
  expect_equal(
    .bincode(x, breaks, right = TRUE, include.lowest = TRUE),
    bin(x, breaks, left_closed = FALSE, include_endpoint = TRUE)
  )
  expect_equal(
    .bincode(x, breaks, right = TRUE, include.lowest = FALSE),
    bin(x, breaks, left_closed = FALSE, include_endpoint = FALSE)
  )

  breaks <- seq(min(x), 5, by = 0.5)

  ### When x is integer here, this might be unexpected result :)

  expect_equal(
    as.integer(.bin(as.integer(x), breaks)),
    bin(as.integer(x), breaks, left_closed = FALSE, codes = FALSE)
  )

  x <- as.double(x)

  expect_equal(
    .bin(x, breaks),
    bin(x, breaks, left_closed = FALSE, codes = FALSE)
  )
  expect_equal(
    .bin(x, breaks, right = FALSE),
    bin(x, breaks, left_closed = TRUE, codes = FALSE)
  )
  expect_equal(
    .bin(x, breaks, include.lowest = TRUE),
    bin(x, breaks, include_endpoint = TRUE, codes = FALSE, left_closed = FALSE)
  )
  expect_equal(
    .bin(x, breaks, right = TRUE, include.lowest = TRUE),
    bin(x, breaks, left_closed = FALSE, include_endpoint = TRUE, codes = FALSE)
  )
  expect_equal(
    .bin(x, breaks, right = TRUE, include.lowest = FALSE),
    bin(x, breaks, left_closed = FALSE, include_endpoint = FALSE, codes = FALSE)
  )

  x <- as.double(-1:10)
  breaks <- 0:11

  expect_equal(
    .bincode(x, breaks, include.lowest = TRUE, right = TRUE),
    bin(x, breaks, include_endpoint = TRUE, left_closed = FALSE)
  )
  expect_equal(
    .bincode(x, breaks, include.lowest = TRUE, right = FALSE),
    bin(x, breaks, include_endpoint = TRUE, left_closed = TRUE)
  )

  expect_equal(
    .bincode(x, c(-Inf, breaks), right = TRUE),
    bin(x, breaks, left_closed = FALSE, include_oob = TRUE)
  )
  expect_equal(
    .bincode(x, c(breaks, Inf), right = FALSE),
    bin(x, breaks, left_closed = TRUE, include_oob = TRUE)
  )

  expect_equal(
    bin(x, breaks, include_oob = TRUE, left_closed = FALSE),
    c(1, 1:11)
  )
  expect_equal(
    bin(x, breaks, include_oob = TRUE, left_closed = TRUE),
    c(NA, 1:11)
  )

  x <- seq(0, 10, 0.5)
  breaks <- seq(1, 9, 0.25)

  expect_equal(
    .bincode(x, breaks),
    bin(x, breaks, left_closed = FALSE)
  )
  expect_equal(
    .bincode(x, breaks, right = FALSE),
    bin(x, breaks, left_closed = TRUE)
  )
  expect_equal(
    .bincode(x, breaks, include.lowest = TRUE),
    bin(x, breaks, include_endpoint = TRUE, left_closed = FALSE)
  )
  expect_equal(
    .bincode(x, breaks, right = TRUE, include.lowest = TRUE),
    bin(x, breaks, left_closed = FALSE, include_endpoint = TRUE)
  )
  expect_equal(
    .bincode(x, breaks, right = TRUE, include.lowest = FALSE),
    bin(x, breaks, left_closed = FALSE, include_endpoint = FALSE)
  )
})

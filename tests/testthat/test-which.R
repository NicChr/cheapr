test_that("c++ which", {
  set.seed(42)
  x <- sample(c(TRUE, FALSE), size = 10^3, replace = TRUE)
  x <- na_insert(x, prop = 0.1)

  expect_identical(which(TRUE), which_(TRUE))
  expect_identical(which(FALSE), which_(FALSE))
  expect_identical(which(logical()), which_(logical()))
  expect_identical(integer(), which_(logical(), invert = TRUE))
  expect_error(which_(1))
  expect_identical(which(x), which_(x))
  expect_identical(which(!x %in% TRUE), which_(x, invert = TRUE))

  x <- c(TRUE, NA, TRUE, NA, NA, rep_len(FALSE, 100))
  expect_identical(which_(x), c(1L, 3L))
  expect_identical(which_(x, invert = TRUE), c(2L, 4L, 5L, 6:length(x)))
})

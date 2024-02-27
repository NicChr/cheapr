test_that("NAs", {
  set.seed(42)
  x <- fill_with_na(1:20, 10)
  x[4] <- 14L
  x[15] <- NA_integer_
  y <- matrix(x, ncol = 4)
  z <- as.data.frame(y)
  names(z) <- paste0("col_", 1:4)
  z

  # Which missing? ----------------------------------------------------------

  expect_identical(which_na(x), which(is.na(x)))
  expect_identical(which_not_na(x), which(!is.na(x)))
  expect_identical(which_na(y), which(is.na(y)))
  expect_identical(which_not_na(y), which(!is.na(y)))
  expect_identical(which_na(z), 5L)
  expect_identical(which_not_na(z), 1:4)


# How many missing? -------------------------------------------------------

  expect_identical(num_na(x), sum(is.na(x)))
  expect_identical(num_na(y), sum(is.na(y)))
  expect_identical(num_na(z), sum(apply(z, 1, function(x) sum(is.na(x)) == ncol(z))))

  expect_error(row_na_counts(x))
  expect_error(col_na_counts(x))
  expect_equal(row_na_counts(y), rowSums(is.na(y)))
  expect_equal(col_na_counts(y), colSums(is.na(y)))
  expect_equal(row_na_counts(z), rowSums(is.na(z)))
  expect_equal(col_na_counts(z), unname(colSums(is.na(z))))


# Any/All NA? -------------------------------------------------------------

  expect_error(row_any_na(x))
  expect_error(col_any_na(x))
  expect_identical(row_any_na(y), apply(y, 1, function(x) sum(is.na(x)) >= 1))
  expect_identical(col_any_na(y), apply(y, 2, function(x) sum(is.na(x)) >= 1))
  expect_identical(row_any_na(z), apply(z, 1, function(x) sum(is.na(x)) >= 1))
  expect_identical(col_any_na(z), unname(apply(z, 2, function(x) sum(is.na(x)) >= 1)))

  expect_error(row_all_na(x))
  expect_error(col_all_na(x))
  expect_identical(row_all_na(y), apply(y, 1, function(x) sum(is.na(x)) >= length(x)))
  expect_identical(col_all_na(y), apply(y, 2, function(x) sum(is.na(x)) >= length(x)))
  expect_identical(row_all_na(z), apply(z, 1, function(x) sum(is.na(x)) >= length(x)))
  expect_identical(col_all_na(z), unname(apply(z, 2, function(x) sum(is.na(x)) >= length(x))))
})

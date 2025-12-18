test_that("whole numbers", {

  expect_identical(is_whole_number(10, na.rm = FALSE), TRUE)
  expect_identical(is_whole_number(10, na.rm = TRUE), TRUE)
  expect_identical(is_whole_number(1.5, na.rm = FALSE), FALSE)
  expect_identical(is_whole_number(1.5, na.rm = TRUE), FALSE)

  expect_identical(is_whole_number(c(NA_real_, 1), na.rm = FALSE), NA)
  expect_identical(is_whole_number(c(NA_real_, 1), na.rm = TRUE), TRUE)


  expect_identical(is_whole_number(c(NA_real_, 1.5), na.rm = FALSE), FALSE)
  expect_identical(is_whole_number(c(NA_real_, 1.5), na.rm = TRUE), FALSE)

  expect_identical(is_whole_number(c(NA_real_, NA_real_), na.rm = FALSE), NA)
  expect_identical(is_whole_number(c(NA_real_, NA_real_), na.rm = TRUE), TRUE)

  # Integers and logical are always true no matter what
  expect_identical(is_whole_number(NA, na.rm = TRUE), TRUE)
  expect_identical(is_whole_number(NA, na.rm = FALSE), TRUE)
  expect_identical(is_whole_number(NA_integer_, na.rm = TRUE), TRUE)
  expect_identical(is_whole_number(NA_integer_, na.rm = FALSE), TRUE)

})

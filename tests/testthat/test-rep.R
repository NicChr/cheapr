test_that("repeating vectors", {

  expect_identical(rep_len_(1:3, 10), rep_len(1:3, 10))
  expect_identical(rep_len_(100, 20), rep_len(100, 20))
  expect_identical(rep_len_(numeric(), 20), rep_len(numeric(), 20))

})

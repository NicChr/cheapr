test_that("reverse", {
  expect_identical(cheapr_rev(numeric()), numeric())
  expect_identical(cheapr_rev(logical()), logical())
  expect_identical(cheapr_rev(character()), character())

  expect_identical(cheapr_rev(1:10), rev(1:10))

  expect_identical(cheapr_rev(sset(iris, 0)), sset(iris, 0))
  expect_identical(cheapr_rev(iris), sset(iris, 150:1))

  expect_identical(cheapr_rev(letters), rev(letters))
})

test_that("reverse", {
  expect_identical(cheapr_rev(numeric()), numeric())
  expect_identical(cheapr_rev(logical()), logical())
  expect_identical(cheapr_rev(character()), character())

  expect_identical(cheapr_rev(1:10), rev(1:10))

  expect_identical(cheapr_rev(sset(iris, 0)), sset(iris, 0))
  expect_identical(cheapr_rev(iris), sset(iris, 150:1))

  expect_identical(cheapr_rev(letters), rev(letters))
})

test_that("local RNG", {

  # Test that with_local_seed always produces the same results
  # given that RNG hasn't been affected between calls

  # Without explicit seed
  expect_identical(
    with_local_seed(rpois(100, 100)),
    with_local_seed(rpois(100, 100))
  )

  # With explicit seed
  expect_identical(
    with_local_seed(rpois(100, 100), .seed = 123),
    with_local_seed(rpois(100, 100), .seed = 123)
  )

  # Test that RNG sequence is unaffected by even nested with_local_seed calls

  set.seed(4712438)
  res <- rnorm(3)
  res2 <- c()
  set.seed(4712438)
  with_local_seed({
    with_local_seed(rnorm(10), .seed = 123)
    rnorm(42)
  }, .seed = 321)
  res2[1] <- rnorm(1)
  with_local_seed({
    with_local_seed(rnorm(32), .seed = 123)
    rnorm(91)
  }, .seed = 321)
  res2[2] <- rnorm(1)
  with_local_seed({
    with_local_seed(rpois(99, 123), .seed = 123)
    runif(77)
  }, .seed = 321)
  res2[3] <- rnorm(1)

  # If these aren't identical then RNG was affected
  # when calling with_local_seed
  expect_identical(res, res2)
})

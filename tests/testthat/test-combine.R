test_that("combining", {

  origin <- as.Date("1970-01-01")

  A <- NULL
  a <- c(FALSE, TRUE, NA)
  b <- 3L
  c <- 7.5
  d <- 1.5 + 3.6i
  e <- c("e", NA_character_)
  f <- origin
  g <- as.POSIXct(origin, tz = "UTC")
  h <- as.POSIXct(origin, tz = "America/New_York")
  i <- factor(c("h", "f", NA), levels = c("h", "H", "f"))
  j <- factor(c("f", "g"), levels = c("f", "G", "g"))
  k <- list(10)
  l <- new_df(x = "ok")

  objs <- list(A, a, b, c, d, e, f, g, h, i, j, k, l)
  result <- Reduce(c_, objs, simplify = FALSE, accumulate = TRUE)

  expect_snapshot(dput(result))
})

test_that("combining", {

  origin <- as.Date("1970-01-01")

  A <- NULL
  a <- c(FALSE, TRUE, NA)
  b <- 3L
  c <- 7.5
  d <- 1.5 + 3.6i
  e <- origin
  f <- as.POSIXct(origin, tz = "UTC")
  g <- as.POSIXct(origin, tz = "America/New_York")
  h <- c("e", NA_character_)
  i <- factor(c("h", "f", NA), levels = c("h", "H", "f"))
  j <- factor(c("f", "g"), levels = c("f", "G", "g"))
  k <- list(10)
  # l <- list(value = 123)
  l <- data.frame(y = "ok")

  objs <- named_list(A, a, b, c, d, e, f, g, h, i, j, k, l)

  result <- vector("list", length(objs))
  names(result) <- names(objs)
  unnamed_objs <- unname(objs)

  expect_equal(foo3(iris), foo3(iris))

  for (ii in seq_along(objs)){
    result <- replace_(
      result, ii, list(a = suppressWarnings(c_(.args = sset(unnamed_objs, seq_len(ii)))))
    )
  }

  expect_snapshot(dput(result))
})

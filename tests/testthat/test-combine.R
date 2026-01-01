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

  expect_true(foo4(1:3))
  expect_true(foo5(1:3))
  expect_equal(foo2(1:3), foo2(1:3))
  expect_equal(alt_class(1:3), alt_class(1:3))
  expect_equal(foobarfoo(1:3), expect_true(bar(1:3)))
  expect_true(bar(1:3))
  expect_true(foobar(1:3))
  expect_true(foo(1:3))
  expect_equal(compact_seq_data(1:3), compact_seq_data(1:3))
  expect_equal(clean_indices(1:3, unnamed_objs, FALSE), clean_indices(1:3, unnamed_objs, FALSE))
  expect_equal(cpp_sset_range(unnamed_objs, 1L, 3L, 1L), cpp_sset_range(unnamed_objs, 1L, 3L, 1L))
  expect_equal(cpp_sset2(unnamed_objs, 1:3, NULL, TRUE, NULL), cpp_sset2(unnamed_objs, 1:3, NULL, TRUE, NULL))
  expect_equal(sset(unnamed_objs, 1:3), sset(unnamed_objs, 1:3))
  expect_equal(sset(unnamed_objs, seq_len(1)), sset(unnamed_objs, seq_len(1)))
  expect_equal(sset(unnamed_objs, 1), sset(unnamed_objs, 1))

  for (ii in seq_along(objs)){
    result <- replace_(
      result, ii, list(a = suppressWarnings(c_(.args = unname(objs)[seq_len(ii)])))
    )
  }

  expect_snapshot(dput(result))
})

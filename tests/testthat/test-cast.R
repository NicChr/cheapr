test_that("casting", {

  origin <- as.Date("1970-01-01")
  a <- 1L
  b <- 2
  c <- FALSE
  d <- 1 + 3i
  e <- "e"
  f <- factor(c("h", "f", NA), levels = c("h", "H", "f"))
  g <- factor(c("f", "g"), levels = c("f", "G", "g"))
  h <- origin
  i <- as.POSIXct(origin, tz = "UTC")
  j <- as.POSIXct(origin, tz = "America/New_York")
  k <- list(10)
  l <- new_df(x = "ok")

  # All permutations of above objects (order matters here)
  all_objs <- new_df(obj = list(a, b, c, d, e, f, g, h, i))

  which_pair <- expand.grid(
    left = seq_len(nrow(all_objs)),
    right = seq_len(nrow(all_objs))
  )

  pairs <- new_df(left = list(), right = list(), .nrows = nrow(which_pair))

  for (i in seq_len(nrow(pairs))){
    pairs <- replace_(pairs, i, new_df(
      left = sset(all_objs, which_pair[[1]][[i]])[[1]],
      right = sset(all_objs, which_pair[[2]][[i]])[[1]]
    ))
  }

  results <- new_list(nrow(pairs))

  for (i in seq_along(results)){
    results[[i]] <- cast_common(pairs[["left"]][[i]], pairs[["right"]][[i]])
  }


  for (i in seq_along(results)){
    results[[i]] <- cast_common(pairs[["left"]][[i]], pairs[["right"]][[i]])
  }

  expect_snapshot(dput(results))
})

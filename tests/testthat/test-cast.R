test_that("casting", {

  origin <- as.Date("1970-01-01")
  a <- 1L
  b <- 2
  c <- FALSE
  d <- 1 + 3i
  e <- "e"
  f <- factor("f")
  g <- as.Date(origin)
  f <- as.POSIXct(g, tz = "Europe/London")
  h <- list(10)
  i <- new_df(x = "ok")

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

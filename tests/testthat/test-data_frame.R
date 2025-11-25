test_that("data frames", {
  expect_identical(new_df(), fast_df())
  expect_identical(new_df(x = 1, y = 2, z = 3), fast_df(x = 1, y = 2, z = 3))
  expect_identical(
    new_df(x = 1, 2, 3, .name_repair = FALSE),
    structure(list(x = 1, 2, 3), row.names = c(NA, -1L), class = "data.frame")
  )
  expect_identical(
    new_df(x = 1, 2, 3),
    structure(list(x = 1, `col_2` = 2, `col_3` = 3), row.names = c(NA, -1L), class = "data.frame")
  )


  # Name repair

  expect_identical(
    new_df(x = 1, y = 2, NULL, x = 3, .name_repair = TRUE),
    list_as_df(list(`x_1` = 1, y = 2, `x_3` = 3))
  )

  # Recycling

  expect_identical(
    new_df(x = 1:3, y = 1:9, .recycle = TRUE),
    new_df(x = rep_len(1:3, 9), y = 1:9)
  )

  # Coercion
  expect_identical(
    as_df(iris),
    iris
  )
  x <- 1:5
  expect_identical(
    as_df(x),
    new_df(value = x)
  )
  expect_identical(
    as_df(matrix(1:10, ncol = 2)),
    new_df(value = 1:10)
  )

  expect_identical(
    as_df(list(x = x)),
    data.frame(x = 1:5)
  )


  # new_df(y = 1, x = matrix(1:10, ncol = 2))
})

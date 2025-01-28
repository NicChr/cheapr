test_that("scalars", {
  a <- with_local_seed(na_insert(round(rnorm(10^3)), prop = 1/3), .seed = 37)

  base_val_rm <- function(x, value){
    x[!x %in% value]
  }

  base_val_replace <- function(x, value, replace){
    x[x %in% value] <- replace
    x
  }

  expect_equal(
    val_count(a, 0), sum(a == 0, na.rm = TRUE)
  )

  expect_equal(
    na_count(a), sum(is.na(a))
  )

  expect_equal(
    val_count(a, NA), na_count(a)
  )

  expect_equal(
   val_rm(a, 0), base_val_rm(a, 0)
  )

  expect_equal(
    val_rm(a, -1), base_val_rm(a, -1)
  )

  expect_equal(
    na_rm(a), a[!is.na(a)]
  )

  expect_equal(
    val_rm(a, NA), na_rm(a)
  )

  expect_equal(
    val_replace(a, -1, 0), base_val_replace(a, -1, 0)
  )
  expect_equal(
    val_replace(a, -1, NA), base_val_replace(a, -1, NA)
  )
  expect_equal(
    val_replace(a, NA, 0), base_val_replace(a, NA, 0)
  )

  expect_equal(
    val_replace(a, NA, 0), na_replace(a, 0)
  )

  expect_equal(val_find(a, NA), which(is.na(a)))
  expect_equal(val_find(a, NA, invert = TRUE), which(!is.na(a)))
  expect_equal(val_find(a, NA, invert = TRUE), na_find(a, invert = TRUE))
  expect_equal(val_find(a, NA), na_find(a))

  expect_equal(na_rm(na_rm(a)), na_rm(a))
  expect_equal(na_rm(a[is.na(a)]), a[0])

  # Data frame empty rows

  expect_identical(na_rm(iris), iris)

  expect_identical(
    na_rm(
      fast_df(
        x = rep(NA, 5),
        y = rep(NA, 5)
      )
    ),
    fast_df(x = logical(), y = logical())
  )

  with_local_seed(
    df <- new_df(x = na_insert(1:20, 7),
                 y = na_insert(1:20, 7))
    , 42
  )

  expect_identical(
    list_as_df(df[rowSums(is.na(df)) < ncol(df), ]),
    na_rm(df)
  )

})

test_that("if else", {


  x <- na_insert(sample.int(3, 100, TRUE), 10)
  z <- complex(real = na_insert(sample(c(0, 1, 2), 100, TRUE), prop = 1/3),
               imaginary = na_insert(sample(c(0, 1, 2), 100, TRUE), prop = 1/3))

  expect_error(cheapr_if_else(x == 3, 1:2, 1))
  expect_error(cheapr_if_else(x == 3, 1, 1:2))
  expect_error(cheapr_if_else(x == 3, 1, 0, na = 1:2))
  expect_error(cheapr_if_else(1L, 1, 0))

  expect_identical(
    cheapr_if_else(x == 3, 0L, 1L),
    ifelse(x == 3, 0L, 1L)
  )
  expect_identical(
    cheapr_if_else(x == 3, 0, 1L),
    ifelse(x == 3, 0, 1L)
  )
  expect_identical(
    cheapr_if_else(x == 3, 0L, 1),
    ifelse(x == 3, 0, 1)
  )

  expect_identical(
    cheapr_if_else(x == 3, .Date(0), .Date(1)),
    .Date(ifelse(x == 3, 0, 1))
  )

  expect_identical(
    cheapr_if_else(x == 3, list(1:3), list(3:1)),
    ifelse(is.na(x), list(NULL), ifelse(x == 3, list(1:3), list(3:1)))
  )

  expect_identical(
    cheapr_if_else(x == 3, list(1:3), list(3:1), na = list(NA_real_)),
    ifelse(is.na(x), list(NA_real_), ifelse(x == 3, list(1:3), list(3:1)))
  )

  expect_identical(
    cheapr_if_else(z == z[1], z[2], z[1]),
    ifelse(is.na(z), z[NA_integer_], ifelse(z == z[1], z[2], z[1]))
  )

  expect_identical(
    cheapr_if_else(x == 3, 0L, "1"),
    ifelse(x == 3, 0, "1")
  )

  expect_identical(
    cheapr_if_else(is.na(x) | !is.na(x), Inf, -Inf),
    ifelse(is.na(x) | !is.na(x), Inf, -Inf)
  )

  expect_identical(
    cheapr_if_else(is.na(x) & !is.na(x), Inf, -Inf),
    ifelse(is.na(x) & !is.na(x), Inf, -Inf)
  )

  expect_identical(
    cheapr_if_else(rep(NA, length(x)), Inf, -Inf),
    as.double(ifelse(rep(NA, length(x)), Inf, -Inf))
  )

  pos_inf <- `class<-`(Inf, "inf_class")
  neg_inf <- `class<-`(-Inf, "inf_class")


  expect_identical(
    cheapr_if_else(is.na(x) | !is.na(x), pos_inf, neg_inf),
    ifelse(is.na(x) | !is.na(x), pos_inf, neg_inf)
  )

  expect_identical(
    cheapr_if_else(is.na(x) & !is.na(x), pos_inf, neg_inf),
    ifelse(is.na(x) & !is.na(x), pos_inf, neg_inf)
  )

  expect_identical(
    cheapr_if_else(rep(NA, length(x)), pos_inf, neg_inf),
    as.double(ifelse(rep(NA, length(x)), pos_inf, neg_inf))
  )

  # expect_identical(
  #   cheapr_if_else(x == 3,
  #                  factor("a", levels = c("a", "c")),
  #                  factor("c", levels = c("c", "b")),
  #                  na = NA)
  # )

  expect_identical(
    cheapr_if_else(is.na(x), iris$Species[1], iris$Species[60]),
    structure(
      as.integer(ifelse(is.na(x), iris$Species[1], iris$Species[60])),
      levels = levels(iris$Species),
      class = "factor"
    )
  )

})

test_that("NAs", {
  set.seed(42)
  x <- fill_with_na(1:20, 10)
  x[4] <- 14L
  x[15] <- NA_integer_
  y <- matrix(x, ncol = 4)
  z <- as.data.frame(y)
  names(z) <- paste0("col_", 1:4)
  z

  expect_identical(num_na(NULL), 0L)
  expect_identical(num_na(numeric()), 0L)

  which_row_any_na <- function(x, invert = FALSE){
    which_(row_any_na(x), invert)
  }
  which_row_all_na <- function(x, invert = FALSE){
    which_(row_all_na(x), invert)
  }
  # Which missing? ----------------------------------------------------------

  expect_identical(which_na(x), which(is.na(x)))
  expect_identical(which_not_na(x), which(!is.na(x)))
  expect_identical(which_na(y), which(is.na(y)))
  expect_identical(which_not_na(y), which(!is.na(y)))
  expect_identical(which_row_all_na(z), 5L)
  expect_identical(which_row_all_na(z, TRUE), 1:4)


# How many missing? -------------------------------------------------------

  expect_identical(num_na(x), sum(is.na(x)))
  expect_identical(num_na(y), sum(is.na(y)))
  expect_identical(num_na(z), sum(is.na(z)))
  # expect_identical(num_na(z), sum(apply(z, 1, function(x) sum(is.na(x)) == ncol(z))))

  expect_error(row_na_counts(x))
  expect_error(col_na_counts(x))
  expect_equal(row_na_counts(y), rowSums(is.na(y)))
  expect_equal(col_na_counts(y), colSums(is.na(y)))
  expect_equal(row_na_counts(z), rowSums(is.na(z)))
  expect_equal(col_na_counts(z), unname(colSums(is.na(z))))


# Any/All NA? -------------------------------------------------------------

  expect_true(any_na(x))
  expect_true(any_na(y))
  expect_true(any_na(z))
  expect_false(all_na(x))
  expect_false(all_na(y))
  expect_false(all_na(z))

  expect_true(any_na(y[5, ]))
  expect_true(all_na(y[5, ]))

  expect_error(row_any_na(x))
  expect_error(col_any_na(x))
  expect_identical(row_any_na(y), apply(y, 1, function(x) sum(is.na(x)) >= 1))
  expect_identical(col_any_na(y), apply(y, 2, function(x) sum(is.na(x)) >= 1))
  expect_identical(row_any_na(z), apply(z, 1, function(x) sum(is.na(x)) >= 1))
  expect_identical(col_any_na(z), unname(apply(z, 2, function(x) sum(is.na(x)) >= 1)))

  expect_error(row_all_na(x))
  expect_error(col_all_na(x))
  expect_identical(row_all_na(y), apply(y, 1, function(x) sum(is.na(x)) >= length(x)))
  expect_identical(col_all_na(y), apply(y, 2, function(x) sum(is.na(x)) >= length(x)))
  expect_identical(row_all_na(z), apply(z, 1, function(x) sum(is.na(x)) >= length(x)))
  expect_identical(col_all_na(z), unname(apply(z, 2, function(x) sum(is.na(x)) >= length(x))))
})

test_that("different classes", {
  set.seed(42)
  a1 <- fill_with_na(1:20, 10)
  a2 <- as.double(a1)
  a3 <- as.character(a1)
  a4 <- complex(real = fill_with_na(1:20, 10),
                imaginary = fill_with_na(1:20, 10))
  l <- list(a1, list(a1, a2), list(a1, a2, a3), list(a1, a2, a3, a4))

  a10 <- a1[0]
  a20 <- a2[0]
  a30 <- a3[0]
  a40 <- a4[0]

  allNA <- function(x){
    all(is.na(x))
  }

  expect_identical(is_na(a10), is.na(a10))
  expect_identical(is_na(a20), is.na(a20))
  expect_identical(is_na(a30), is.na(a30))
  expect_identical(is_na(a40), is.na(a40))

  expect_identical(is_na(a1), is.na(a1))
  expect_identical(is_na(a2), is.na(a2))
  expect_identical(is_na(a3), is.na(a3))
  expect_identical(is_na(a4), is.na(a4))

  expect_identical(num_na(a10), sum(is.na(a10)))
  expect_identical(num_na(a20), sum(is.na(a20)))
  expect_identical(num_na(a30), sum(is.na(a30)))
  expect_identical(num_na(a40), sum(is.na(a40)))

  expect_identical(num_na(a1), sum(is.na(a1)))
  expect_identical(num_na(a2), sum(is.na(a2)))
  expect_identical(num_na(a3), sum(is.na(a3)))
  expect_identical(num_na(a4), sum(is.na(a4)))

  expect_identical(num_na(l[0]), sum(is.na(unlist(l[0]))))
  expect_identical(num_na(l), sum(is.na(unlist(l))))

  expect_identical(which_na(a10), which(is.na(a10)))
  expect_identical(which_na(a20), which(is.na(a20)))
  expect_identical(which_na(a30), which(is.na(a30)))
  expect_identical(which_na(a40), which(is.na(a40)))

  expect_identical(which_na(a1), which(is.na(a1)))
  expect_identical(which_na(a2), which(is.na(a2)))
  expect_identical(which_na(a3), which(is.na(a3)))
  expect_identical(which_na(a4), which(is.na(a4)))

  expect_identical(which_not_na(a10), which(!is.na(a10)))
  expect_identical(which_not_na(a20), which(!is.na(a20)))
  expect_identical(which_not_na(a30), which(!is.na(a30)))
  expect_identical(which_not_na(a40), which(!is.na(a40)))

  expect_identical(which_not_na(a1), which(!is.na(a1)))
  expect_identical(which_not_na(a2), which(!is.na(a2)))
  expect_identical(which_not_na(a3), which(!is.na(a3)))
  expect_identical(which_not_na(a4), which(!is.na(a4)))

  expect_identical(any_na(a10), anyNA(a10))
  expect_identical(any_na(a20), anyNA(a20))
  expect_identical(any_na(a30), anyNA(a30))
  expect_identical(any_na(a40), anyNA(a40))

  expect_identical(any_na(a1), anyNA(a1))
  expect_identical(any_na(a2), anyNA(a2))
  expect_identical(any_na(a3), anyNA(a3))
  expect_identical(any_na(a4), anyNA(a4))

  expect_identical(all_na(a10), allNA(a10))
  expect_identical(all_na(a20), allNA(a20))
  expect_identical(all_na(a30), allNA(a30))
  expect_identical(all_na(a40), allNA(a40))

  expect_identical(all_na(a1), allNA(a1))
  expect_identical(all_na(a2), allNA(a2))
  expect_identical(all_na(a3), allNA(a3))
  expect_identical(all_na(a4), allNA(a4))
})

test_that("multiple cores", {
  options(cheapr.cores = 2)
  set.seed(42)
  a1 <- fill_with_na(1:10^5, 10)
  a2 <- as.double(a1)
  a3 <- as.character(a1)
  a4 <- complex(real = fill_with_na(1:10^5, 10),
                imaginary = fill_with_na(1:10^5, 10))
  l <- list(a1, list(a1, a2), list(a1, a2, a3), list(a1, a2, a3, a4))

  allNA <- function(x){
    all(is.na(x))
  }
  expect_identical(is_na(a1), is.na(a1))
  expect_identical(is_na(a2), is.na(a2))
  expect_identical(is_na(a3), is.na(a3))
  expect_identical(is_na(a4), is.na(a4))

  expect_identical(num_na(a1), sum(is.na(a1)))
  expect_identical(num_na(a2), sum(is.na(a2)))
  expect_identical(num_na(a3), sum(is.na(a3)))
  expect_identical(num_na(a4), sum(is.na(a4)))

  expect_identical(num_na(l), sum(is.na(unlist(l))))

  expect_identical(which_na(a1), which(is.na(a1)))
  expect_identical(which_na(a2), which(is.na(a2)))
  expect_identical(which_na(a3), which(is.na(a3)))
  expect_identical(which_na(a4), which(is.na(a4)))

  expect_identical(which_not_na(a1), which(!is.na(a1)))
  expect_identical(which_not_na(a2), which(!is.na(a2)))
  expect_identical(which_not_na(a3), which(!is.na(a3)))
  expect_identical(which_not_na(a4), which(!is.na(a4)))

  expect_identical(any_na(a1), anyNA(a1))
  expect_identical(any_na(a2), anyNA(a2))
  expect_identical(any_na(a3), anyNA(a3))
  expect_identical(any_na(a4), anyNA(a4))

  expect_identical(all_na(a1), allNA(a1))
  expect_identical(all_na(a2), allNA(a2))
  expect_identical(all_na(a3), allNA(a3))
  expect_identical(all_na(a4), allNA(a4))
  options(cheapr.cores = 1)
})

test_that("lists", {

  # NA counts of empty row/col data frames
  expect_identical(
    row_na_counts(data.frame(row.names = 1:100)),
    integer(100)
  )
  expect_identical(
    row_na_counts(data.frame(x = 1:100,
                             y = rnorm(100))[0, , drop = FALSE]),
    integer()
  )
  expect_identical(
    col_na_counts(data.frame(row.names = 1:100)),
    integer()
  )
  expect_identical(
    col_na_counts(data.frame(x = 1:100,
                             y = rnorm(100))[0, , drop = FALSE]),
    integer(2)
  )
  # any NA empty row/col data frames
  expect_identical(
    row_any_na(data.frame(row.names = 1:100)),
    logical(100)
  )
  expect_identical(
    row_any_na(data.frame(x = rep(NA, 100),
                             y = rnorm(100))[0, , drop = FALSE]),
    logical(0)
  )
  expect_identical(
    col_any_na(data.frame(row.names = 1:100)),
    logical(0)
  )
  expect_identical(
    col_any_na(data.frame(x = rep(NA, 100),
                             y = rep(NA, 100))[0, , drop = FALSE]),
    logical(2)
  )

  # all NA empty row/col data frames
  expect_identical(
    row_all_na(data.frame(row.names = 1:100)),
    logical(100)
  )
  expect_identical(
    row_all_na(data.frame(x = rep(NA, 100),
                          y = rep(NA, 100))[0, , drop = FALSE]),
    logical(0)
  )
  expect_identical(
    col_all_na(data.frame(row.names = 1:100)),
    logical(0)
  )
  expect_identical(
    col_all_na(data.frame(x = rep(NA, 100),
                          y = rep(NA, 100))[0, , drop = FALSE]),
    logical(2)
  )

  x <- list(1, 1:5, NA, list(1:10, list(c(2, NA, NA))))
  expect_identical(num_na(x), 3L)
  expect_true(any_na(x))
  expect_false(all_na(x))

  df <- list(x = list(1, 1:5, NA, list(NA, NA), integer()), y = rep(NA, 5))
  attributes(df) <- list(class = "data.frame",
                         row.names = c(NA_integer_, -5),
                         names = as.character(names(df)))
  df
  expect_identical(
    row_na_counts(df),
    vapply(df$x, function(x) (sum(is.na(unlist(x))) > 0L) + 1L, 0L)
  )
  expect_identical(
    col_na_counts(df),
    c(2L, 5L)
  )

  ## Triggering parallel code
  options(cheapr.cores = 2)
  set.seed(912389)
  df <- list(a = as.list(fill_with_na(sample(letters, 10^5 + 1, TRUE), n = 10^4)),
             b = fill_with_na(rnorm(10^5 + 1), 10^3),
             c = fill_with_na(sample.int(100, 10^5 + 1, TRUE), 10^3),
             d = complex(real = fill_with_na(sample.int(100, 10^5 + 1, TRUE), 10^3),
                         imaginary = fill_with_na(sample.int(100, 10^5 + 1, TRUE), 10^3)),
             e = fill_with_na(sample(letters, 10^5 + 1, TRUE), n = 10^4))

  attributes(df) <- list(class = "data.frame",
                         row.names = c(NA_integer_, -as.integer(10^5 + 1)),
                         names = as.character(names(df)))
  expect_identical(
    row_na_counts(df),
    # Add 1 every time ALL unlist elements are NA
    vapply(df$a, function(x) (sum(is.na(unlist(x))) > 0L), 0L) +
      is.na(df$b) +
      is.na(df$c) +
      is.na(df$d) +
      is.na(df$e)
  )
  expect_identical(
    col_na_counts(df),
    c(10000L, 1000L, 1000L, sum(is.na(df$d)), 10000L)
  )
  options(cheapr.cores = 1)
})

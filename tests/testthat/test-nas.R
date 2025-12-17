
is.na2 <- function(x){
  is.na(x)
}

recursive_na_count <- function(x){
  sum(unlist(rapply(x, \(x) sum(is.na2(x)), how = "list")))
}

all_na2 <- function(x){
  all(is.na2(x))
}
any_na2 <- function(x){
  any(is.na2(x))
}

test_that("NAs", {
  set.seed(42)
  x <- na_insert(1:20, 10)
  x[4] <- 14L
  x[15] <- NA_integer_
  y <- matrix(x, ncol = 4)
  z <- as.data.frame(y)
  names(z) <- paste0("col_", 1:4)
  z

  expect_identical(na_count(NULL), 0L)
  expect_identical(na_count(numeric()), 0L)

  which_row_any_na <- function(x, invert = FALSE){
    which_(row_any_na(x), invert)
  }
  which_row_all_na <- function(x, invert = FALSE){
    which_(row_all_na(x), invert)
  }
  # Which missing? ----------------------------------------------------------

  expect_identical(which_na(x), which(is.na2(x)))
  expect_identical(which_not_na(x), which(!is.na2(x)))
  expect_identical(which_na(y), which(is.na2(y)))
  expect_identical(which_not_na(y), which(!is.na2(y)))
  expect_identical(which_row_all_na(z), 5L)
  expect_identical(which_row_all_na(z, TRUE), 1:4)


# How many missing? -------------------------------------------------------

  expect_identical(na_count(x), sum(is.na2(x)))
  expect_identical(na_count(y), sum(is.na2(y)))
  expect_identical(na_count(z), sum(is.na2(z)))
  # expect_identical(na_count(z), sum(apply(z, 1, function(x) sum(is.na2(x)) == ncol(z))))

  expect_error(row_na_counts(x))
  expect_error(col_na_counts(x))
  expect_equal(row_na_counts(y), rowSums(is.na2(y)))
  expect_equal(col_na_counts(y), colSums(is.na2(y)))
  expect_equal(row_na_counts(z), rowSums(is.na2(z)))
  expect_equal(col_na_counts(z), unname(colSums(is.na2(z))))


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
  expect_identical(row_any_na(y), apply(y, 1, function(x) sum(is.na2(x)) >= 1))
  expect_identical(col_any_na(y), apply(y, 2, function(x) sum(is.na2(x)) >= 1))
  expect_identical(row_any_na(z), apply(z, 1, function(x) sum(is.na2(x)) >= 1))
  expect_identical(col_any_na(z), unname(apply(z, 2, function(x) sum(is.na2(x)) >= 1)))

  expect_error(row_all_na(x))
  expect_error(col_all_na(x))
  expect_identical(row_all_na(y), apply(y, 1, function(x) sum(is.na2(x)) >= length(x)))
  expect_identical(col_all_na(y), apply(y, 2, function(x) sum(is.na2(x)) >= length(x)))
  expect_identical(row_all_na(z), apply(z, 1, function(x) sum(is.na2(x)) >= length(x)))
  expect_identical(col_all_na(z), unname(apply(z, 2, function(x) sum(is.na2(x)) >= length(x))))
})

test_that("different classes", {
  set.seed(42)
  a1 <- na_insert(1:20, 10)
  a2 <- as.double(a1)
  a3 <- as.character(a1)
  a4 <- complex(real = na_insert(1:20, 10),
                imaginary = na_insert(1:20, 10))
  l <- list(a1, list(a1, a2), list(a1, a2, a3), list(a1, a2, a3, a4))

  a10 <- a1[0]
  a20 <- a2[0]
  a30 <- a3[0]
  a40 <- a4[0]

  expect_identical(is_na(a10), is.na2(a10))
  expect_identical(is_na(a20), is.na2(a20))
  expect_identical(is_na(a30), is.na2(a30))
  expect_identical(is_na(a40), is.na2(a40))

  expect_identical(is_na(a1), is.na2(a1))
  expect_identical(is_na(a2), is.na2(a2))
  expect_identical(is_na(a3), is.na2(a3))
  expect_identical(is_na(a4), is.na2(a4))

  expect_identical(na_count(a10), sum(is.na2(a10)))
  expect_identical(na_count(a20), sum(is.na2(a20)))
  expect_identical(na_count(a30), sum(is.na2(a30)))
  expect_identical(na_count(a40), sum(is.na2(a40)))

  expect_identical(na_count(a1), sum(is.na2(a1)))
  expect_identical(na_count(a2), sum(is.na2(a2)))
  expect_identical(na_count(a3), sum(is.na2(a3)))
  expect_identical(na_count(a4), sum(is.na2(a4)))

  expect_identical(na_count(l[0]), recursive_na_count(l[0]))
  expect_identical(na_count(l), recursive_na_count(l))

  expect_identical(which_na(a10), which(is.na2(a10)))
  expect_identical(which_na(a20), which(is.na2(a20)))
  expect_identical(which_na(a30), which(is.na2(a30)))
  expect_identical(which_na(a40), which(is.na2(a40)))

  expect_identical(which_na(a1), which(is.na2(a1)))
  expect_identical(which_na(a2), which(is.na2(a2)))
  expect_identical(which_na(a3), which(is.na2(a3)))
  expect_identical(which_na(a4), which(is.na2(a4)))

  expect_identical(which_not_na(a10), which(!is.na2(a10)))
  expect_identical(which_not_na(a20), which(!is.na2(a20)))
  expect_identical(which_not_na(a30), which(!is.na2(a30)))
  expect_identical(which_not_na(a40), which(!is.na2(a40)))

  expect_identical(which_not_na(a1), which(!is.na2(a1)))
  expect_identical(which_not_na(a2), which(!is.na2(a2)))
  expect_identical(which_not_na(a3), which(!is.na2(a3)))
  expect_identical(which_not_na(a4), which(!is.na2(a4)))

  expect_identical(any_na(a10), any_na2(a10))
  expect_identical(any_na(a20), any_na2(a20))
  expect_identical(any_na(a30), any_na2(a30))
  expect_identical(any_na(a40), any_na2(a40))

  expect_identical(any_na(a1), any_na2(a1))
  expect_identical(any_na(a2), any_na2(a2))
  expect_identical(any_na(a3), any_na2(a3))
  expect_identical(any_na(a4), any_na2(a4))

  expect_identical(all_na(a10), all_na2(a10))
  expect_identical(all_na(a20), all_na2(a20))
  expect_identical(all_na(a30), all_na2(a30))
  expect_identical(all_na(a40), all_na2(a40))

  expect_identical(all_na(a1), all_na2(a1))
  expect_identical(all_na(a2), all_na2(a2))
  expect_identical(all_na(a3), all_na2(a3))
  expect_identical(all_na(a4), all_na2(a4))
})

test_that("multiple cores", {
  set.seed(42)
  a1 <- na_insert(1:10^5, 10)
  a2 <- as.double(a1)
  a3 <- as.character(a1)
  a4 <- complex(real = na_insert(1:10^5, 10),
                imaginary = na_insert(1:10^5, 10))
  l <- list(a1, list(a1, a2), list(a1, a2, a3), list(a1, a2, a3, a4))

  expect_identical(is_na(a1), is.na2(a1))
  expect_identical(is_na(a2), is.na2(a2))
  expect_identical(is_na(a3), is.na2(a3))
  expect_identical(is_na(a4), is.na2(a4))

  expect_identical(na_count(a1), sum(is.na2(a1)))
  expect_identical(na_count(a2), sum(is.na2(a2)))
  expect_identical(na_count(a3), sum(is.na2(a3)))
  expect_identical(na_count(a4), sum(is.na2(a4)))

  expect_identical(na_count(l), recursive_na_count(l))

  expect_identical(which_na(a1), which(is.na2(a1)))
  expect_identical(which_na(a2), which(is.na2(a2)))
  expect_identical(which_na(a3), which(is.na2(a3)))
  expect_identical(which_na(a4), which(is.na2(a4)))

  expect_identical(which_not_na(a1), which(!is.na2(a1)))
  expect_identical(which_not_na(a2), which(!is.na2(a2)))
  expect_identical(which_not_na(a3), which(!is.na2(a3)))
  expect_identical(which_not_na(a4), which(!is.na2(a4)))

  expect_identical(any_na(a1), any_na2(a1))
  expect_identical(any_na(a2), any_na2(a2))
  expect_identical(any_na(a3), any_na2(a3))
  expect_identical(any_na(a4), any_na2(a4))

  expect_identical(all_na(a1), all_na2(a1))
  expect_identical(all_na(a2), all_na2(a2))
  expect_identical(all_na(a3), all_na2(a3))
  expect_identical(all_na(a4), all_na2(a4))
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
    rep(TRUE, 100)
    # logical(100)
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
    rep(TRUE, 2)
    # logical(2)
  )

  x <- list(1, 1:5, NA, list(1:10, list(c(2, NA, NA))))
  expect_identical(na_count(x), 3L)
  expect_true(any_na(x))
  expect_false(all_na(x))

  df <- list(x = list(1, 1:5, NA, list(NA, NA), integer()), y = rep(NA, 5))
  attributes(df) <- list(class = "data.frame",
                         row.names = c(NA_integer_, -5),
                         names = as.character(names(df)))
  df
  expect_identical(
    row_na_counts(df),
    vapply(df$x, function(x) (sum(is.na2(unlist(x))) > 0L) + 1L, 0L)
  )
  expect_identical(
    col_na_counts(df),
    c(2L, 5L)
  )

  set.seed(912389)
  df <- list(a = as.list(na_insert(sample(letters, 10^5 + 1, TRUE), n = 10^4)),
             b = na_insert(rnorm(10^5 + 1), 10^3),
             c = na_insert(sample.int(100, 10^5 + 1, TRUE), 10^3),
             d = complex(real = na_insert(sample.int(100, 10^5 + 1, TRUE), 10^3),
                         imaginary = na_insert(sample.int(100, 10^5 + 1, TRUE), 10^3)),
             e = na_insert(sample(letters, 10^5 + 1, TRUE), n = 10^4))

  attributes(df) <- list(class = "data.frame",
                         row.names = c(NA_integer_, -as.integer(10^5 + 1)),
                         names = as.character(names(df)))
  expect_identical(
    row_na_counts(df),
    # Add 1 every time ALL unlist elements are NA
    vapply(df$a, is_na, TRUE) +
      is.na2(df$b) +
      is.na2(df$c) +
      is.na2(df$d) +
      is.na2(df$e)
  )
  expect_identical(
    col_na_counts(df),
    c(10000L, 1000L, 1000L, sum(is.na2(df$d)), 10000L)
  )
})

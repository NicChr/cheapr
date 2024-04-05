options(cheapr.cores = 2)
test_that("subsetting", {
  set.seed(37) # Totally random number :)
  a <- rnorm(10^3)
  a[sample.int(10^3, 10)] <- Inf
  a[sample.int(10^3, 10)] <- -Inf
  a[sample.int(10^3, 10)] <- NaN
  b <- sample(-100:100, 10^3, TRUE)
  c <- sample(letters, 10^3, TRUE)
  d <- complex(real = rnorm(10^3),
               imaginary = rnorm(10^3))
  e <- .Date(b)
  f <- .POSIXct(b)
  g <- as.POSIXlt(f)
  h <- as.factor(c)
  # i <- as.list(a)
  i <- vector("list", 10^3)
  for (k in seq_along(i)){
    i[[k]] <- c(a[k], b[k])
  }
  i
  j <- vapply(sample(letters, 10^3, TRUE), charToRaw, raw(1))
  k <- a > 0

  i1 <- sample.int(10^3)
  i2 <- sample.int(100, 10^4, TRUE)
  i3 <- i2
  i3[sample.int(length(i3), 100)] <- 0L
  i3[sample.int(length(i3), 100)] <- NA
  i4 <- -sample.int(100, 10^4, TRUE)
  i4[sample.int(length(i4), 100)] <- 0L
  i5 <- sample(0:5e03, 10^4, TRUE)
  i6 <- 0L
  i7 <- NA_integer_
  i8 <- sample.int(10^3, 10^5, TRUE)
  i9 <- 3:10^3
  i10 <- 1e03:3
  i11 <- 1:10^3
  i12 <- 11:500
  i12 <- 1:333
  i13 <- 2e03:11
  i14 <- 1e03:1
  i15 <- 450:1
  i16 <- -1:-200
  i17 <- -200:-1
  i18 <- 0:-233
  i19 <- -233:0
  i20 <- -255:-555
  i21 <- -500:-3000
  i22 <- -2000:-3000
  i23 <- 2000:3000
  i24 <- -1:-1000
  i25 <- -111:-1000

  # i7 <- NA # This doesn't match
  # i8 <- NA_integer_

  objs_to_test <- letters[1:10]
  ind_to_test <- paste0("i", 1:25)

  for (obj in objs_to_test){
    if (!is.raw(get(obj))){
      assign(obj, fill_with_na(get(obj), n = 111))
    }
  }

  df <- data.frame(a, b, c, d, e, f, g, h)
  df$i <- i
  df$j <- j
  df$k <- k

  test_df <- expand.grid(objs_to_test, ind_to_test, stringsAsFactors = FALSE)
  names(test_df) <- c("obj", "ind")
  # test_df <- test_df |> dplyr::filter(obj != "g")
  for (i in seq_len(nrow(test_df))){
    r_obj <- get(test_df$obj[i])
    r_ind <- get(test_df$ind[i])
    expect_identical(
     sset(r_obj, r_ind),
     r_obj[r_ind]
    )
  }
  # Base data frame subset except row names are reset
  base_sset <- function(...){
    out <- `[`(...)
    attr(out, "row.names") <- .set_row_names(nrow(out))
    out
  }
  for (ind in ind_to_test){
    r_ind <- get(ind)
    expect_identical(
     sset(df, r_ind),
     base_sset(df, r_ind, 1:ncol(df), drop = FALSE)
     # vctrs::vec_slice(df, r_ind)
    )
  }

  empty_df <- df[0, , drop = FALSE]

  expect_identical(
    lapply(empty_df, `[`, 0:10),
    lapply(empty_df, sset, 0:10)
  )
  expect_identical(
    lapply(df, `[`, 0),
    lapply(df, sset, 0)
  )
  expect_identical(
    lapply(empty_df, `[`, NA_integer_),
    lapply(empty_df, sset, NA_integer_)
  )

  expect_identical(
    sset(empty_df, 0:20),
    base_sset(empty_df, 0:20, 1:ncol(empty_df), drop = FALSE)
  )
})

test_that("errors", {
  expect_null(sset(NULL))
  expect_null(sset(NULL, 1:10))
  expect_error(sset(iris$Sepal.Length, c(-5, 5)))
  expect_error(sset(iris$Sepal.Length, 10:-10))
  expect_error(sset(globalenv()))
})

test_that("misc", {
  expect_identical(
    cpp_sset_range(10:1, 3, 0, -1),
    rev(cpp_sset_range(10:1, 0, 3, 1))
  )
  expect_identical(cpp_sset_range(1:10, 0, 0, 1), integer())
  expect_identical(cpp_sset_range(integer(), 0, 0, 1), integer())
  expect_identical(cpp_sset_range(letters, 0, 0, 1), character())
  expect_identical(cpp_sset_range(as.list(letters), 0, 0, 1), list())
  expect_identical(cpp_sset_range(as.list(letters)[0], 0, 0, -1), list())
  expect_identical(cpp_sset_range(letters[0], 0, 0, -1), character())
  expect_error(cpp_sset_range(letters, 1, 10, 2))
  expect_error(cpp_sset_range(letters, 1, 10, -1))
  expect_error(cpp_sset_range(global(), 1, 10, 1))
})

# test_that("fatal error", {
#   set.seed(43)
#   i <- as.integer(sample(seq(0, 1000, 100), 10^6, TRUE))
#   x <- flights$day
#   # x <- sample.int(30, 100000, TRUE)
#   cpp_sset(x, -i)
# }
# )

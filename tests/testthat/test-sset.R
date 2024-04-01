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
  # i7 <- NA # This doesn't match
  # i8 <- NA_integer_

  objs_to_test <- letters[1:8]
  ind_to_test <- paste0("i", 1:8)

  for (obj in objs_to_test){
    assign(obj, fill_with_na(get(obj), n = 111))
  }

  df <- data.frame(a, b, c, d, e, f, g, h)

  test_df <- expand.grid(objs_to_test, ind_to_test, stringsAsFactors = FALSE)
  names(test_df) <- c("obj", "ind")


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

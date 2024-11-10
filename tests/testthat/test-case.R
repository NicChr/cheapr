test_that("matching", {

  with_local_seed({


    x <- floor(rnorm(10^4, 100, 50))
    x <- na_insert(x, prop = 1/5)

    keys <- c("zero", "one", "two", "three", "four", "five", "six", "seven",
              "eight", "nine", "ten", "NA", "Other")
    values <- c(0:10, NA)

    target <- keys[match(x, values, nomatch = length(values) + 1L)]
    result <- val_match(
      x, 0 ~ "zero",
      1 ~ "one",
      2 ~ "two",
      3 ~ "three",
      4 ~ "four",
      5 ~ "five",
      6 ~ "six",
      7 ~ "seven",
      8 ~ "eight",
      9 ~ "nine",
      10 ~ "ten",
      NA ~ "NA",
      .default = "Other"
    )

    expect_identical(result, target)
    expect_identical(
      cheapr_table(target),
      c(Other = 7870L, `NA` = 2000L, eight = 6L, seven = 19L, nine = 19L,
        one = 10L, zero = 10L, six = 11L, five = 11L, two = 10L, ten = 8L,
        four = 12L, three = 14L)
    )


    expect_identical(
      case(x == 0 ~ "zero",
        x == 1 ~ "one",
        x == 2 ~ "two",
        x == 3 ~ "three",
        x == 4 ~ "four",
        x == 5 ~ "five",
        x == 6 ~ "six",
        x == 7 ~ "seven",
        x == 8 ~ "eight",
        x == 9 ~ "nine",
        x == 10 ~ "ten",
        is.na(x) ~ "NA",
        TRUE ~ "Other"
      ),
      target
    )

  },
  .seed = 12345)
})


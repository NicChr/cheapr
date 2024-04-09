test_that("sequences", {
  expect_error(seq_(1, 0, by = 1))
  expect_error(seq_size(1, 0, by = 1))
  expect_identical(sequence_(numeric()), integer())
  expect_identical(sequence(1:3), sequence_(1:3))
  expect_identical(c(1, 1, 1.1, 1, 1.1, 1.2),
                   sequence_(1:3, by = 0.1))

  expect_error(sequence_(1, integer(), integer()))
  expect_error(sequence_(1, numeric(), numeric()))
  expect_identical(sequence_(integer(), integer(), integer()), integer())
  expect_identical(sequence_(numeric(), numeric(), numeric()), numeric())

  # Numeric vs integer
  expect_equal(sequence(c(3, 2), by = c(-1, 1)),
               sequence_(c(3, 2), by = c(-1, 1)))

  set.seed(98376234)
  x <- sample(0:99)
  from <- sample(11:20)
  by <- sample(1:5)
  expect_identical(sequence(x, from, by), sequence_(x, from, by))
  expect_equal(sequence(x, from, by), sequence_(x, as.double(from), as.double(by)))

  expect_identical(sequence(123, from = 1L, by = 0L),
                   sequence_(123, from = 1L, by = 0L))

  expect_identical(sequence_(0, from = 1L, by = 0L),
                   integer())
  expect_identical(sequence_(integer(), from = 1L, by = 0L),
                   integer())

  expect_equal(
    seq_(from = Sys.Date() + c(0, 11),
         to = Sys.Date() + c(0, 12, 13, 23),
         by = c(2, 7)),
    sequence(seq_size(from = Sys.Date() + c(0, 11),
                      to = Sys.Date() + c(0, 12, 13, 23),
                      by = c(2, 7)),
             from = Sys.Date() + c(0, 11),
             by = c(2, 7))
  )

  expect_identical(seq_(-1L, 10L, 1L), seq.int(-1L, 10L, 1L))
  expect_equal(seq_(0, 20, 0.2), seq.int(0, 20, 0.2))

  expect_identical(
    window_sequence(c(5, 0, 1, 4), 3, add_id = TRUE),
    c(`1` = 1L, `1` = 2L, `1` = 3L, `1` = 3L, `1` = 3L, `3` = 1L,
      `4` = 1L, `4` = 2L, `4` = 3L, `4` = 3L)
  )
  expect_identical(
    lag_sequence(c(5, 0, 1, 4), 3, add_id = TRUE),
    c(`1` = 0L, `1` = 1L, `1` = 2L, `1` = 3L, `1` = 3L, `3` = 0L,
      `4` = 0L, `4` = 1L, `4` = 2L, `4` = 3L),
  )
  expect_identical(
    lag_sequence(c(5, 0, 1, 4), 3, add_id = TRUE, partial = FALSE),
    c(`1` = NA, `1` = NA, `1` = NA, `1` = 3L, `1` = 3L, `3` = NA,
      `4` = NA, `4` = NA, `4` = NA, `4` = 3L),
  )
  expect_identical(
    lead_sequence(c(5, 0, 1, 4), 3, add_id = TRUE),
    c(`1` = 3L, `1` = 3L, `1` = 2L, `1` = 1L, `1` = 0L, `3` = 0L,
      `4` = 3L, `4` = 2L, `4` = 1L, `4` = 0L),
  )
  expect_identical(
    lead_sequence(c(5, 0, 1, 4), 3, add_id = TRUE, partial = FALSE),
    c(`1` = 3L, `1` = 3L, `1` = NA, `1` = NA, `1` = NA, `3` = NA,
      `4` = 3L, `4` = NA, `4` = NA, `4` = NA),
  )

})


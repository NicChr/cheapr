test_that("math operations", {
  options(cheapr.cores = 2)

  # The first 2 will not trigger multi-threaded calculations
  make_test_data1 <- function(){
    set.seed(3742)
    assign("x", fill_with_na(rnorm(10^3), 50),
           envir = parent.frame())
  }
  make_test_data2 <- function(){
    set.seed(3742)
    assign("x", fill_with_na(sample.int(100, 10^3, TRUE), 50),
           envir = parent.frame())
  }
  make_test_data3 <- function(){
    set.seed(3742)
    assign("x", fill_with_na(rnorm(10^5), 10^3),
           envir = parent.frame())
  }
  make_test_data4 <- function(){
    set.seed(3742)
    assign("x", fill_with_na(sample.int(100, 10^5, TRUE), 10^3),
           envir = parent.frame())
  }

  expect_error(set_abs(iris))
  expect_error(set_floor(iris))
  expect_error(set_ceiling(iris))
  expect_error(set_trunc(iris))
  expect_error(set_log(iris))
  expect_error(set_exp(iris))
  expect_error(set_sqrt(iris))
  expect_error(set_round(iris))

  make_test_data1()
  expect_equal(
   round(x, 2),
   set_round(sset(x), 2)
  )
  make_test_data1()
  expect_equal(
    round(x),
    set_round(sset(x))
  )
  make_test_data1()
  expect_equal(round(x, 1), set_round(sset(x), 1))
  make_test_data1()
  expect_equal(abs(x), set_abs(sset(x)))
  make_test_data1()
  expect_identical(floor(x), set_floor(sset(x)))
  make_test_data1()
  expect_identical(ceiling(x), set_ceiling(sset(x)))
  make_test_data1()
  expect_identical(trunc(x), set_trunc(sset(x)))
  make_test_data1()
  expect_equal(exp(x), set_exp(sset(x)))
  make_test_data1()
  expect_equal(log(abs(x)), set_log(set_abs(sset(x))))
  make_test_data1()
  expect_equal(log10(abs(x)), set_log(set_abs(sset(x)), base = 10))
  make_test_data1()
  expect_equal(sqrt(abs(x)), set_sqrt(set_abs(sset(x))))
  make_test_data1()
  expect_equal(-x, set_change_sign(sset(x)))


  make_test_data2()
  expect_equal(round(x, 1),set_round(x, 1))
  make_test_data2()
  expect_equal(abs(x), set_abs(sset(x)))
  make_test_data2()
  expect_equal(floor(x), set_floor(sset(x)))
  make_test_data2()
  expect_equal(ceiling(x), set_ceiling(sset(x)))
  make_test_data2()
  expect_equal(trunc(x), set_trunc(sset(x)))
  make_test_data2()
  expect_equal(exp(x), suppressWarnings(set_exp(sset(x))))
  make_test_data2()
  expect_equal(log(abs(x)), suppressWarnings(set_log(set_abs(sset(x)))))
  make_test_data2()
  expect_equal(log10(abs(x)), suppressWarnings(set_log(set_abs(sset(x)), base = 10)))
  make_test_data2()
  expect_equal(sqrt(abs(x)), suppressWarnings(set_sqrt(set_abs(sset(x)))))
  make_test_data2()
  expect_equal(-x, set_change_sign(sset(x)))

  make_test_data3()
  expect_equal(round(x, 1), set_round(sset(x), 1))
  make_test_data3()
  expect_equal(abs(x), set_abs(sset(x)))
  make_test_data3()
  expect_equal(floor(x), set_floor(sset(x)))
  make_test_data3()
  expect_equal(ceiling(x), set_ceiling(sset(x)))
  make_test_data3()
  expect_equal(trunc(x), set_trunc(sset(x)))
  make_test_data3()
  expect_equal(exp(x), set_exp(sset(x)))
  make_test_data3()
  expect_equal(log(abs(x)), set_log(set_abs(sset(x))))
  make_test_data3()
  expect_equal(log10(abs(x)), set_log(set_abs(sset(x)), base = 10))
  make_test_data3()
  expect_equal(sqrt(abs(x)), set_sqrt(set_abs(sset(x))))
  make_test_data3()
  expect_equal(-x, set_change_sign(sset(x)))

  make_test_data4()
  expect_equal(round(x, 1),set_round(sset(x), 1))
  make_test_data4()
  expect_equal(abs(x), set_abs(sset(x)))
  make_test_data4()
  expect_equal(floor(x), set_floor(sset(x)))
  make_test_data4()
  expect_equal(ceiling(x), set_ceiling(sset(x)))
  make_test_data4()
  expect_equal(trunc(x), set_trunc(sset(x)))
  make_test_data4()
  expect_equal(exp(x), suppressWarnings(set_exp(sset(x))))
  make_test_data4()
  expect_equal(log(abs(x)), suppressWarnings(set_log(set_abs(sset(x)))))
  make_test_data4()
  expect_equal(log10(abs(x)), suppressWarnings(set_log(set_abs(sset(x)), base = 10)))
  make_test_data4()
  expect_equal(sqrt(abs(x)), suppressWarnings(set_sqrt(set_abs(sset(x)))))
  make_test_data4()
  expect_equal(-x, set_change_sign(sset(x)))
  options(cheapr.cores = 1)
})

test_that("more math operations", {
  options(cheapr.cores = 2)

  # The first 2 will not trigger multi-threaded calculations
  make_test_data1 <- function(){
    set.seed(3742)
    assign("x", fill_with_na(rnorm(10^3), 50),
           envir = parent.frame())
    assign("y", c(0L, NA_integer_, 3:10), envir = parent.frame())
    assign("z", sample(c(Inf, -Inf, NA_real_, sequence_(100 - 3, 0, 0.1))),
           envir = parent.frame())
  }
  make_test_data2 <- function(){
    set.seed(3742)
    assign("x", fill_with_na(sample.int(100, 10^3, TRUE), 50),
           envir = parent.frame())
    assign("y", c(0L, NA_integer_, 3:10), envir = parent.frame())
    assign("z", sample(c(Inf, -Inf, NA_real_, sequence_(100 - 3, 0, 0.1))),
           envir = parent.frame())
  }
  make_test_data3 <- function(){
    set.seed(3742)
    assign("x", fill_with_na(rnorm(10^5), 10^3),
           envir = parent.frame())
    assign("y", c(0L, NA_integer_, 3:10), envir = parent.frame())
    assign("z", sample(c(Inf, -Inf, NA_real_, sequence_(1000 - 3, 0, 0.1))),
           envir = parent.frame())
  }
  make_test_data4 <- function(){
    set.seed(3742)
    assign("x", fill_with_na(sample.int(100, 10^5, TRUE), 10^3),
           envir = parent.frame())
    assign("y", c(0L, NA_integer_, 3:10), envir = parent.frame())
    assign("z", sample(c(Inf, -Inf, NA_real_, sequence_(1000 - 3, 0, 0.1))),
           envir = parent.frame())
  }

  expect_error(set_add(iris))
  expect_error(set_subtract(iris))
  expect_error(set_multiply(iris))
  expect_error(set_divide(iris))

  make_test_data1()
  expect_equal(
    x + y,
    set_add(sset(x), y)
  )
  expect_equal(
    x + z,
    set_add(sset(x), z)
  )
  expect_equal(
    x - y,
    set_subtract(sset(x), y)
  )
  expect_equal(
    x - z,
    set_subtract(sset(x), z)
  )
  expect_equal(
    x * y,
    set_multiply(sset(x), y)
  )
  expect_equal(
    x * z,
    set_multiply(sset(x), z)
  )
  expect_equal(
    x / y,
    set_divide(sset(x), y)
  )
  expect_equal(
    x / z,
    set_divide(sset(x), z)
  )
  expect_equal(
    x^y,
    set_pow(sset(x), y)
  )
  expect_equal(
    round(x, y),
    set_round(sset(x), y)
  )
  # expect_equal(
  #   round(x, as.integer(z)),
  #   set_round(sset(x), z)
  # )
  expect_equal(
    suppressWarnings(log(x, base = y)),
    set_log(sset(x), y)
  )
  expect_equal(
    suppressWarnings(log(x, base = z)),
    set_log(sset(x), z)
  )

  make_test_data2()
  expect_equal(
    x + y,
    set_add(sset(x), y)
  )
  expect_equal(
    x + z,
    suppressWarnings(set_add(sset(x), z))
  )
  expect_equal(
    x - y,
    set_subtract(sset(x), y)
  )
  expect_equal(
    x - z,
    suppressWarnings(set_subtract(sset(x), z))
  )
  expect_equal(
    x * y,
    set_multiply(sset(x), y)
  )
  expect_equal(
    x * z,
    suppressWarnings(set_multiply(sset(x), z))
  )
  expect_equal(
    x / y,
    suppressWarnings(set_divide(sset(x), y))
  )
  expect_equal(
    x / z,
    suppressWarnings(set_divide(sset(x), z))
  )
  expect_equal(
    x^y,
    suppressWarnings(set_pow(sset(x), y))
  )
  expect_equal(
    x,
    # round(x, y),
    set_round(sset(x), y)
  )
  expect_equal(
    suppressWarnings(log(x, base = y)),
    suppressWarnings(set_log(sset(x), y))
  )

  make_test_data3()
  expect_equal(
    x + y,
    set_add(sset(x), y)
  )
  expect_equal(
    x + z,
    set_add(sset(x), z)
  )
  expect_equal(
    x - y,
    set_subtract(sset(x), y)
  )
  expect_equal(
    x - z,
    set_subtract(sset(x), z)
  )
  expect_equal(
    x * y,
    set_multiply(sset(x), y)
  )
  expect_equal(
    x * z,
    set_multiply(sset(x), z)
  )
  expect_equal(
    x / y,
    set_divide(sset(x), y)
  )
  expect_equal(
    x^y,
    set_pow(sset(x), y)
  )
  expect_equal(
    round(x, y),
    set_round(sset(x), y)
  )
  expect_equal(
    suppressWarnings(log(x, base = y)),
    suppressWarnings(set_log(sset(x), y))
  )

  make_test_data4()
  expect_equal(
    x + y,
    set_add(sset(x), y)
  )
  expect_equal(
    x + z,
    suppressWarnings(set_add(sset(x), z))
  )
  expect_equal(
    x - y,
    set_subtract(sset(x), y)
  )
  expect_equal(
    x - z,
    suppressWarnings(set_subtract(sset(x), z))
  )
  expect_equal(
    x * y,
    set_multiply(sset(x), y)
  )
  expect_equal(
    x * z,
    suppressWarnings(set_multiply(sset(x), z))
  )
  expect_equal(
    x / y,
    suppressWarnings(set_divide(sset(x), y))
  )
  expect_equal(
    x / z,
    suppressWarnings(set_divide(sset(x), z))
  )
  expect_equal(
    x^y,
    suppressWarnings(set_pow(sset(x), y))
  )
  expect_equal(
    # round(x, y),
    x,
    set_round(sset(x), y)
  )
  expect_equal(
    suppressWarnings(log(x, base = y)),
    suppressWarnings(set_log(sset(x), y))
  )
  options(cheapr.cores = 1)
})

test_that("math operations", {
  options(cheapr.cores = 2)

  # The first 2 will not trigger multi-threaded calculations
  make_test_data1 <- function(){
    set.seed(3742)
    assign("x", na_insert(rnorm(10^3), 50),
           envir = parent.frame())
  }
  make_test_data2 <- function(){
    set.seed(3742)
    assign("x", na_insert(sample.int(100, 10^3, TRUE), 50),
           envir = parent.frame())
  }
  make_test_data3 <- function(){
    set.seed(3742)
    assign("x", na_insert(rnorm(10^5), 10^3),
           envir = parent.frame())
  }
  make_test_data4 <- function(){
    set.seed(3742)
    assign("x", na_insert(sample.int(100, 10^5, TRUE), 10^3),
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
   base::round(x, 2),
   set_round(r_copy(x), 2)
  )
  make_test_data1()
  expect_equal(
    base::round(x),
    set_round(r_copy(x))
  )
  make_test_data1()
  expect_equal(base::round(x, 1), set_round(r_copy(x), 1))

  expect_equal(abs(x), set_abs(r_copy(x)))

  expect_identical(floor(x), set_floor(r_copy(x)))

  expect_identical(ceiling(x), set_ceiling(r_copy(x)))

  expect_identical(base::trunc(x), set_trunc(r_copy(x)))

  expect_equal(exp(x), set_exp(r_copy(x)))

  expect_equal(log(abs(x)), set_log(set_abs(r_copy(x))))

  expect_equal(log10(abs(x)), set_log(set_abs(r_copy(x)), base = 10))

  expect_equal(sqrt(abs(x)), set_sqrt(set_abs(r_copy(x))))

  expect_equal(-x, set_change_sign(r_copy(x)))


  make_test_data2()
  expect_equal(base::round(x, 1), set_round(x, 1))

  expect_equal(abs(x), set_abs(r_copy(x)))

  expect_equal(floor(x), set_floor(r_copy(x)))

  expect_equal(ceiling(x), set_ceiling(r_copy(x)))

  expect_equal(base::trunc(x), set_trunc(r_copy(x)))

  expect_equal(exp(x), suppressWarnings(set_exp(r_copy(x))))

  expect_equal(log(abs(x)), suppressWarnings(set_log(set_abs(r_copy(x)))))

  expect_equal(log10(abs(x)), suppressWarnings(set_log(set_abs(r_copy(x)), base = 10)))

  expect_equal(sqrt(abs(x)), suppressWarnings(set_sqrt(set_abs(r_copy(x)))))

  expect_equal(-x, set_change_sign(r_copy(x)))

  make_test_data3()
  expect_equal(base::round(x, 1), set_round(r_copy(x), 1))

  expect_equal(abs(x), set_abs(r_copy(x)))

  expect_equal(floor(x), set_floor(r_copy(x)))

  expect_equal(ceiling(x), set_ceiling(r_copy(x)))

  expect_equal(base::trunc(x), set_trunc(r_copy(x)))

  expect_equal(exp(x), set_exp(r_copy(x)))

  expect_equal(log(abs(x)), set_log(set_abs(r_copy(x))))

  expect_equal(log10(abs(x)), set_log(set_abs(r_copy(x)), base = 10))

  expect_equal(sqrt(abs(x)), set_sqrt(set_abs(r_copy(x))))

  expect_equal(-x, set_change_sign(r_copy(x)))

  make_test_data4()
  expect_equal(base::round(x, 1),set_round(r_copy(x), 1))

  expect_equal(abs(x), set_abs(r_copy(x)))

  expect_equal(floor(x), set_floor(r_copy(x)))

  expect_equal(ceiling(x), set_ceiling(r_copy(x)))

  expect_equal(base::trunc(x), set_trunc(r_copy(x)))

  expect_equal(exp(x), suppressWarnings(set_exp(r_copy(x))))

  expect_equal(log(abs(x)), suppressWarnings(set_log(set_abs(r_copy(x)))))

  expect_equal(log10(abs(x)), suppressWarnings(set_log(set_abs(r_copy(x)), base = 10)))

  expect_equal(sqrt(abs(x)), suppressWarnings(set_sqrt(set_abs(r_copy(x)))))

  expect_equal(-x, set_change_sign(r_copy(x)))
  options(cheapr.cores = 1)
})

test_that("more math operations", {
  options(cheapr.cores = 2)

  # The first 2 will not trigger multi-threaded calculations
  make_test_data1 <- function(){
    set.seed(3742)
    assign("x", na_insert(rnorm(10^3), 50),
           envir = parent.frame())
    assign("y", c(0L, NA_integer_, 3:10), envir = parent.frame())
    assign("z", sample(c(Inf, -Inf, NA_real_, sequence_(100 - 3, 0, 0.1))),
           envir = parent.frame())
  }
  make_test_data2 <- function(){
    set.seed(3742)
    assign("x", na_insert(sample.int(100, 10^3, TRUE), 50),
           envir = parent.frame())
    assign("y", c(0L, NA_integer_, 3:10), envir = parent.frame())
    assign("z", sample(c(Inf, -Inf, NA_real_, sequence_(100 - 3, 0, 0.1))),
           envir = parent.frame())
  }
  make_test_data3 <- function(){
    set.seed(3742)
    assign("x", na_insert(rnorm(10^5), 10^3),
           envir = parent.frame())
    assign("y", c(0L, NA_integer_, 3:10), envir = parent.frame())
    assign("z", sample(c(Inf, -Inf, NA_real_, sequence_(1000 - 3, 0, 0.1))),
           envir = parent.frame())
  }
  make_test_data4 <- function(){
    set.seed(3742)
    assign("x", na_insert(sample.int(100, 10^5, TRUE), 10^3),
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
    set_add(r_copy(x), y)
  )
  expect_equal(
    x + z,
    set_add(r_copy(x), z)
  )
  expect_equal(
    x - y,
    set_subtract(r_copy(x), y)
  )
  expect_equal(
    x - z,
    set_subtract(r_copy(x), z)
  )
  expect_equal(
    x * y,
    set_multiply(r_copy(x), y)
  )
  expect_equal(
    x * z,
    set_multiply(r_copy(x), z)
  )
  expect_equal(
    x / y,
    set_divide(r_copy(x), y)
  )
  expect_equal(
    x / z,
    set_divide(r_copy(x), z)
  )
  expect_equal(
    x^y,
    set_pow(r_copy(x), y)
  )
  expect_equal(
    base::round(x, y),
    set_round(r_copy(x), y)
  )
  # expect_equal(
  #   base::round(x, as.integer(z)),
  #   set_round(r_copy(x), z)
  # )
  expect_equal(
    suppressWarnings(log(x, base = y)),
    set_log(r_copy(x), y)
  )
  expect_equal(
    suppressWarnings(log(x, base = z)),
    set_log(r_copy(x), z)
  )

  make_test_data2()
  expect_equal(
    x + y,
    set_add(r_copy(x), y)
  )
  expect_equal(
    x + z,
    suppressWarnings(set_add(r_copy(x), z))
  )
  expect_equal(
    x - y,
    set_subtract(r_copy(x), y)
  )
  expect_equal(
    x - z,
    suppressWarnings(set_subtract(r_copy(x), z))
  )
  expect_equal(
    x * y,
    set_multiply(r_copy(x), y)
  )
  expect_equal(
    x * z,
    suppressWarnings(set_multiply(r_copy(x), z))
  )
  expect_equal(
    x / y,
    suppressWarnings(set_divide(r_copy(x), y))
  )
  expect_equal(
    x / z,
    suppressWarnings(set_divide(r_copy(x), z))
  )
  expect_equal(
    x^y,
    suppressWarnings(set_pow(r_copy(x), y))
  )
  expect_equal(
    x,
    # base::round(x, y),
    set_round(r_copy(x), y)
  )
  expect_equal(
    suppressWarnings(log(x, base = y)),
    suppressWarnings(set_log(r_copy(x), y))
  )

  make_test_data3()
  expect_equal(
    x + y,
    set_add(r_copy(x), y)
  )
  expect_equal(
    x + z,
    set_add(r_copy(x), z)
  )
  expect_equal(
    x - y,
    set_subtract(r_copy(x), y)
  )
  expect_equal(
    x - z,
    set_subtract(r_copy(x), z)
  )
  expect_equal(
    x * y,
    set_multiply(r_copy(x), y)
  )
  expect_equal(
    x * z,
    set_multiply(r_copy(x), z)
  )
  expect_equal(
    x / y,
    set_divide(r_copy(x), y)
  )
  expect_equal(
    x^y,
    set_pow(r_copy(x), y)
  )
  expect_equal(
    base::round(x, y),
    set_round(r_copy(x), y)
  )
  expect_equal(
    suppressWarnings(log(x, base = y)),
    suppressWarnings(set_log(r_copy(x), y))
  )

  make_test_data4()
  expect_equal(
    x + y,
    set_add(r_copy(x), y)
  )
  expect_equal(
    x + z,
    suppressWarnings(set_add(r_copy(x), z))
  )
  expect_equal(
    x - y,
    set_subtract(r_copy(x), y)
  )
  expect_equal(
    x - z,
    suppressWarnings(set_subtract(r_copy(x), z))
  )
  expect_equal(
    x * y,
    set_multiply(r_copy(x), y)
  )
  expect_equal(
    x * z,
    suppressWarnings(set_multiply(r_copy(x), z))
  )
  expect_equal(
    x / y,
    suppressWarnings(set_divide(r_copy(x), y))
  )
  expect_equal(
    x / z,
    suppressWarnings(set_divide(r_copy(x), z))
  )
  expect_equal(
    x^y,
    suppressWarnings(set_pow(r_copy(x), y))
  )
  expect_equal(
    # base::round(x, y),
    x,
    set_round(r_copy(x), y)
  )
  expect_equal(
    suppressWarnings(log(x, base = y)),
    suppressWarnings(set_log(r_copy(x), y))
  )
  options(cheapr.cores = 1)
})

test_that("zero-length vectors", {
  x <- 1
  y <- numeric()

  expect_error(set_add(x, y))
  expect_error(set_subtract(x, y))
  expect_error(set_divide(x, y))
  expect_error(set_multiply(x, y))
  expect_error(set_log(x, y))
  expect_error(set_pow(x, y))
  expect_error(set_round(x, y))

  expect_identical(set_add(y, x), numeric())
  expect_identical(set_subtract(y, x), numeric())
  expect_identical(set_divide(y, x), numeric())
  expect_identical(set_multiply(y, x), numeric())
  expect_identical(set_log(y, x), numeric())
  expect_identical(set_pow(y, x), numeric())
  expect_identical(set_round(y, x), numeric())

  x <- 1L
  y <- integer()

  expect_error(set_add(x, y))
  expect_error(set_subtract(x, y))
  expect_error(set_divide(x, y))
  expect_error(set_multiply(x, y))
  expect_error(set_log(x, y))
  expect_error(set_pow(x, y))
  expect_error(set_round(x, y))

  expect_identical(set_add(y, x), integer())
  expect_identical(set_subtract(y, x), integer())
  expect_identical(set_divide(as.double(y), x), numeric())
  expect_identical(set_multiply(y, x), integer())
  expect_identical(set_log(as.double(y), x), numeric())
  expect_identical(set_pow(as.double(y), x), numeric())
  expect_identical(set_round(y, x), integer())
})

test_that("math operations", {

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

  expect_error(abs_(iris))
  expect_error(floor_(iris))
  expect_error(ceiling_(iris))
  expect_error(trunc_(iris))
  expect_error(log_(iris))
  expect_error(exp_(iris))
  expect_error(sqrt_(iris))
  expect_error(round_(iris))

  make_test_data1()
  expect_equal(
    base::round(x, 2),
    round_(x, 2)
  )
  make_test_data1()
  expect_equal(
    base::round(x),
    round_(x)
  )
  make_test_data1()
  expect_equal(base::round(x, 1), round_(x, 1))

  expect_equal(abs(x), abs_(x))

  expect_identical(floor(x), floor_(x))

  expect_identical(ceiling(x), ceiling_(x))

  expect_identical(base::trunc(x), trunc_(x))

  expect_equal(exp(x), exp_(x))

  expect_equal(log(abs(x)), log_(abs_(x)))

  expect_equal(log10(abs(x)), log_(abs_(x), base = 10))

  expect_equal(sqrt(abs(x)), sqrt_(abs_(x)))

  expect_equal(-x, negate_(x))


  make_test_data2()
  expect_equal(base::round(x, 1), round_(x, 1))

  expect_equal(abs(x), abs_(x))

  expect_equal(floor(x), floor_(x))

  expect_equal(ceiling(x), ceiling_(x))

  expect_equal(base::trunc(x), trunc_(x))

  expect_equal(exp(x), exp_(x))

  expect_equal(log(abs(x)), log_(abs_(x)))

  expect_equal(log10(abs(x)), log_(abs_(x), base = 10))

  expect_equal(sqrt(abs(x)), sqrt_(abs_(x)))

  expect_equal(-x, negate_(x))

  make_test_data3()
  expect_equal(base::round(x, 1), round_(x, 1))

  expect_equal(abs(x), abs_(x))

  expect_equal(floor(x), floor_(x))

  expect_equal(ceiling(x), ceiling_(x))

  expect_equal(base::trunc(x), trunc_(x))

  expect_equal(exp(x), exp_(x))

  expect_equal(log(abs(x)), log_(abs_(x)))

  expect_equal(log10(abs(x)), log_(abs_(x), base = 10))

  expect_equal(sqrt(abs(x)), sqrt_(abs_(x)))

  expect_equal(-x, negate_(x))

  make_test_data4()
  expect_equal(base::round(x, 1),round_(x, 1))

  expect_equal(abs(x), abs_(x))

  expect_equal(floor(x), floor_(x))

  expect_equal(ceiling(x), ceiling_(x))

  expect_equal(base::trunc(x), trunc_(x))

  expect_equal(exp(x), exp_(x))

  expect_equal(log(abs(x)), log_(abs_(x)))

  expect_equal(log10(abs(x)), log_(abs_(x), base = 10))

  expect_equal(sqrt(abs(x)), sqrt_(abs_(x)))

  expect_equal(-x, negate_(x))
})

test_that("more math operations", {

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

  expect_error(add_(iris))
  expect_error(subtract_(iris))
  expect_error(multiply_(iris))
  expect_error(divide_(iris))

  make_test_data1()
  expect_equal(
    x + y,
    add_(x, y)
  )
  expect_equal(
    x + z,
    add_(x, z)
  )
  expect_equal(
    x - y,
    subtract_(x, y)
  )
  expect_equal(
    x - z,
    subtract_(x, z)
  )
  expect_equal(
    x * y,
    multiply_(x, y)
  )
  expect_equal(
    x * z,
    multiply_(x, z)
  )
  expect_equal(
    x / y,
    divide_(x, y)
  )
  expect_equal(
    x / z,
    divide_(x, z)
  )
  expect_equal(
    x^y,
    pow_(x, y)
  )
  expect_equal(
    base::round(x, y),
    round_(x, y)
  )
  expect_equal(
    suppressWarnings(log(x, base = y)),
    log_(x, y)
  )
  expect_equal(
    suppressWarnings(log(x, base = z)),
    log_(x, z)
  )

  make_test_data2()
  expect_equal(
    x + y,
    add_(x, y)
  )
  expect_equal(
    x + z,
    add_(x, z)
  )
  expect_equal(
    x - y,
    subtract_(x, y)
  )
  expect_equal(
    x - z,
    subtract_(x, z)
  )
  expect_equal(
    x * y,
    multiply_(x, y)
  )
  expect_equal(
    x * z,
    multiply_(x, z)
  )
  expect_equal(
    x / y,
    divide_(x, y)
  )
  expect_equal(
    x / z,
    divide_(x, z)
  )
  expect_equal(
    x^y,
    pow_(x, y)
  )
  expect_equal(
    x,
    # base::round(x, y),
    round_(x, y)
  )
  expect_equal(
    suppressWarnings(log(x, base = y)),
    log_(x, y)
  )

  make_test_data3()
  expect_equal(
    x + y,
    add_(x, y)
  )
  expect_equal(
    x + z,
    add_(x, z)
  )
  expect_equal(
    x - y,
    subtract_(x, y)
  )
  expect_equal(
    x - z,
    subtract_(x, z)
  )
  expect_equal(
    x * y,
    multiply_(x, y)
  )
  expect_equal(
    x * z,
    multiply_(x, z)
  )
  expect_equal(
    x / y,
    divide_(x, y)
  )
  expect_equal(
    x^y,
    pow_(x, y)
  )
  expect_equal(
    base::round(x, y),
    round_(x, y)
  )
  expect_equal(
    suppressWarnings(log(x, base = y)),
    suppressWarnings(log_(x, y))
  )

  make_test_data4()
  expect_equal(
    x + y,
    add_(x, y)
  )
  expect_equal(
    x + z,
    add_(x, z)
  )
  expect_equal(
    x - y,
    subtract_(x, y)
  )
  expect_equal(
    x - z,
    subtract_(x, z)
  )
  expect_equal(
    x * y,
    multiply_(x, y)
  )
  expect_equal(
    x * z,
    multiply_(x, z)
  )
  expect_equal(
    x / y,
    divide_(x, y)
  )
  expect_equal(
    x / z,
    divide_(x, z)
  )
  expect_equal(
    x^y,
    pow_(x, y)
  )
  expect_equal(
    # base::round(x, y),
    x,
    round_(x, y)
  )
  expect_equal(
    suppressWarnings(log(x, base = y)),
    log_(x, y)
  )
})

test_that("zero-length vectors", {
  x <- 1
  y <- numeric()

  expect_identical(add_(x, y), numeric())
  expect_identical(subtract_(x, y), numeric())
  expect_identical(divide_(x, y), numeric())
  expect_identical(multiply_(x, y), numeric())
  expect_identical(log_(x, y), numeric())
  expect_identical(pow_(x, y), numeric())
  expect_identical(round_(x, y), numeric())

  expect_identical(add_(y, x), numeric())
  expect_identical(subtract_(y, x), numeric())
  expect_identical(divide_(y, x), numeric())
  expect_identical(multiply_(y, x), numeric())
  expect_identical(log_(y, x), numeric())
  expect_identical(pow_(y, x), numeric())
  expect_identical(round_(y, x), numeric())

  x <- 1L
  y <- integer()

  expect_identical(add_(x, y), integer())
  expect_identical(subtract_(x, y), integer())
  expect_identical(divide_(as.double(x), y), numeric())
  expect_identical(multiply_(x, y), integer())
  expect_identical(log_(as.double(x), y), numeric())
  expect_identical(pow_(as.double(x), y), numeric())
  expect_identical(round_(x, y), integer())

  expect_identical(add_(y, x), integer())
  expect_identical(subtract_(y, x), integer())
  expect_identical(divide_(as.double(y), x), numeric())
  expect_identical(multiply_(y, x), integer())
  expect_identical(log_(as.double(y), x), numeric())
  expect_identical(pow_(as.double(y), x), numeric())
  expect_identical(round_(y, x), integer())
})

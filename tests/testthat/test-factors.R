test_that("factors", {
  expect_identical(factor_(NULL), factor())

  x <- rnorm(10)
  y <- na_insert(x, 5)
  expect_identical(factor(x), factor_(x))
  expect_identical(factor(y), factor_(y))
  expect_identical(
    factor(y, exclude = NULL),
    factor_(y, na_exclude = FALSE)
  )

  airquality <- datasets::airquality

  # Levels sorted by order of first appearance
  expect_identical(
    factor_(airquality$Ozone, order = FALSE),
    factor(airquality$Ozone, levels = unique(airquality$Ozone))
  )
  expect_identical(
    factor_(airquality$Ozone, order = FALSE, na_exclude = FALSE),
    factor(airquality$Ozone, levels = unique(airquality$Ozone),
           exclude = NULL)
  )

  # Explicit levels
  expect_identical(
    factor(airquality$Day, levels = 10:15),
    factor_(airquality$Day, levels = 10:15)
  )
  expect_identical(
    factor(airquality$Day, levels = c(NA, 10:15)),
    factor_(airquality$Day, levels = c(NA, 10:15))
  )
  # Used and unused levels
  fct <- factor_(airquality$Ozone, levels = c(NA, 10:100), na_exclude = FALSE)
  expect_identical(
    levels_used(fct),
    intersect(levels(fct), airquality$Ozone)
  )
  expect_identical(
    levels_unused(fct),
    setdiff(levels(fct), airquality$Ozone)
  )
  expect_identical(
    levels_factor(factor_(airquality$Day)),
    as.factor(as.integer(levels(factor_(airquality$Day))))
  )

  # Datetimes
  now <- structure(1709021560.32021, class = c("POSIXct", "POSIXt"))
  x <- seq(now, now + 8800 * 60 * 60, by = 3600)
  expect_identical(factor_(x[1:10]), factor(x[1:10]))
})

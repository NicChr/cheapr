test_that("overview", {
  options(cheapr.digits = 2)
  expect_snapshot(overview(airquality, hist = FALSE))
  iris2 <- iris
  iris2$large <- iris2$Sepal.Length >= 6
  iris2$Species2 <- as.character(iris2$Species)
  iris2 <- iris2[which(iris2$Species %in% unique(iris2$Species)[1:2]), , drop = FALSE]
  expect_snapshot(overview(iris, hist = FALSE))
  expect_snapshot(overview(iris2, hist = FALSE))
  expect_snapshot(overview(warpbreaks, hist = FALSE))
  expect_snapshot(overview(ToothGrowth, hist = FALSE))
  expect_identical(overview(ts(1:10))$numeric,
                   overview(data.frame(x = ts(1:10)))$time_series)
  # expect_snapshot(overview(EuStockMarkets, hist = TRUE))
})

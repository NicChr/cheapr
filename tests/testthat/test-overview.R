test_that("overview", {
  expect_snapshot(overview(airquality, hist = TRUE))
  iris2 <- iris
  iris2$large <- iris2$Sepal.Length >= 6
  iris2$Species2 <- as.character(iris2$Species)
  iris2 <- iris2[which(iris2$Species %in% unique(iris2$Species)[1:2]), , drop = FALSE]
  expect_snapshot(overview(iris, hist = TRUE))
  expect_snapshot(overview(iris2, hist = TRUE))
  expect_snapshot(overview(warpbreaks, hist = TRUE))
  expect_snapshot(overview(ToothGrowth, hist = TRUE))
})

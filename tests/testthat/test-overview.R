test_that("overview", {
  options(cheapr.digits = 2)
  expect_snapshot(overview(numeric(), hist = TRUE))
  expect_snapshot(overview(numeric(), hist = FALSE))
  expect_snapshot(overview(airquality, hist = FALSE))
  iris2 <- iris
  iris2$large <- iris2$Sepal.Length >= 6
  iris2$Species2 <- as.character(iris2$Species)
  iris2 <- iris2[which(iris2$Species %in% unique(iris2$Species)[1:2]), , drop = FALSE]
  expect_snapshot(overview(iris, hist = FALSE))
  expect_snapshot(overview(iris2, hist = FALSE))
  expect_snapshot(overview(warpbreaks, hist = FALSE))
  expect_snapshot(overview(ToothGrowth, hist = FALSE))
  expect_identical(overview(ts(1:10)),
                   overview(data.frame(x = ts(1:10))))
  x <- c(1.14258786651076, -0.364521829529459, 0.657994899969532, -1.23613869620319,
         0.817144011037811, 1.32624440032353, 1.29713129390224, -0.389266962953551,
         0.80756330472006, 0.767873979361665, -0.105195294912731, 0.643003584192831,
         -2.52978196144816, 0.243452492963461, -1.59622496369583, 0.0758326451931188,
         1.59003458564609, 0.52747123961581, -0.699964248520901, -0.44330920205848,
         -1.23317865134766, -2.3631591543418, 0.858336607111522, 1.50765939182219,
         0.0275182238317008)
  df <- data.frame(
    y = ts(x),
    x = ts(x)
  ) |>
    col_c(
      as.data.frame(matrix(rep(x, 5), ncol = 5)) |>
        lapply(ts) |>
        new_df(.args = _) |>
        stats::setNames(paste0("z_Series ", 1:5))
    )
  expect_snapshot(overview(df))
  expect_snapshot(overview(ts(matrix(x, ncol = 5))))
  expect_snapshot(overview(EuStockMarkets))
})

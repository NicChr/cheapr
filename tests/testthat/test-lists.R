test_that("list assignment", {
  l <- as.list(iris)

  expect_equal(
    list_assign(l, `names<-`(new_list(length(l)), names(l))),
    `names<-`(list(), character())
  )

  expect_equal(
    list_assign(l, list(x = 1L, Species = NULL, y = 0L, y = NULL)),
    c(l[c(1, 2, 3, 4)], list(x = 1L, y = 0L))
  )

})

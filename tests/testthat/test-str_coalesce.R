test_that("normal tests", {
  # Normal examples
  expect_equal(str_coalesce("", "hello"), "hello")
  expect_equal(str_coalesce("", NA, "goodbye"), "goodbye")

  # '' always preferred
  expect_equal(str_coalesce("", NA), "")
  expect_equal(str_coalesce(NA, ""), "")

  # Unless there are only NAs
  expect_identical(str_coalesce(NA, NA), NA_character_)

  # `str_coalesce` is vectorised

  with_local_seed({
    x <- val_insert(letters, "", n = 10)
    y <- val_insert(LETTERS, "", n = 10)
  }, 42)

  expect_equal(str_coalesce(x, y), ifelse(x != "" & !is.na(x), x, y))

  expect_equal(
    str_coalesce(NULL, letters, letters, NULL, NULL, letters, NULL, NULL),
    letters
  )
})

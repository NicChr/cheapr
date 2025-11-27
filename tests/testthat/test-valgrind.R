test_that("memory leak test", {
  expect_error(cheapr_do_memory_leak_test())
})

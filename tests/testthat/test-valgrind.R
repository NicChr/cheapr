test_that("memory leak test", {
  expect_error(cheapr_do_memory_leak_test()) # No memory leak should happen here
  expect_error(cheapr_unsafe_init_memory_leak()) # This SHOULD cause a memory leak
})

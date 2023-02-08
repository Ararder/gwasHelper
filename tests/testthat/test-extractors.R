test_that("The environmental variable gwasHelper_repo is set", {
  expect_true(Sys.getenv("gwasHelper_repo") != "")
})

test_that("The environmental variable gwasHelper_repo is set", {
  expect_true(fs::dir_exists(Sys.getenv("gwasHelper_repo")))
})

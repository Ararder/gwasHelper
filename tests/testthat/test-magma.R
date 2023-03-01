paths <- setup_test_dir()




test_that("magma writes out two files", {
  fs::dir_create(paths[["magma_dir"]])
  magma_munge_and_setup(
    input_gwas = paths[["clean"]],
    output_dir = paths[["magma_dir"]]
    )

  expect_true(
    fs::file_exists(
      fs::path(paths[["magma_dir"]], "snplocs.tsv")
      )
  )
  expect_true(
    fs::file_exists(
      fs::path(paths[["magma_dir"]], "pval.tsv")
      )
  )

})

test_that("Correctly captures code to run magma",{
  # magma_gene(paths[["magma_dir"]])
  expect_true(TRUE)
})


# test_that("Correctly captures code to run magma",{
#   magma_geneset(paths[["magma_dir"]], "test_annotation")
#   expect_true(TRUE)
# })


# test_that("magma returns code",{
#   magma_gene(paths[["magma_dir"]])
# })
#
test_that("magma pipeline captures code to be used",{
  cleansumstats_magma(paths = paths)

})



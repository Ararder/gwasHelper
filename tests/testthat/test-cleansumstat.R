# setups a temp dir that mimics the structure of gwasHelper_repo
testdir <- withr::local_tempdir()
withr::local_envvar("gwasHelper_repo" = testdir)
paths <- filepath_manager(test_path("data/gwas_repo/bipolar_disorder.txt.gz"), "bipolar2021")
fs::dir_create(paths[["base"]])
col_map <- c(
  "col_SNP: SNP",
  "col_POS: BP",
  "col_BETA: Effect",
  "col_SE: SE",
  "col_P: Pval",
  "col_N: 'N'"
)
clean_sumstats(
  infile = test_path("data/gwas_repo/bipolar_disorder.txt.gz"),
  col_map = col_map,
  name = "bipolar2021"
)



test_that("clean_sumstats_job writes a meta.txt file to the base directory directory", {
  make_cleansumstats_job(paths, model = "logistic", col_map = col_map)
  expect_true(
    yaml::read_yaml(paths[["metafile"]])["col_SNP"] == "SNP"
  )
})

test_that("clean_sumstats writes out a slurm script", {
  expect_true(
    fs::file_exists(paths[["base_job"]])
  )
})

test_that("clean_sumstats_job creates the directory for pldsc stiletti", {
  expect_true(
   fs::dir_exists(paths[["pldsc_siletti"]])
  )
})





setup_test_dir <- function() {
  testdir <- withr::local_tempdir(clean=FALSE)
  withr::local_envvar("gwasHelper_repo" = testdir)
  paths <- filepath_manager(test_path("data/gwas_repo/bipolar_disorder.txt.gz"), "bipolar2021")
  fs::dir_create(paths[["base"]])
  fs::dir_create(paths[["cleaned_dir"]])


  fs::file_copy(test_path("data/cleaned_GRCh38.gz"), paths[["clean"]])
  fs::file_copy(test_path("data/cleaned_GRCh37.gz"), paths[["build37"]])

  withr::defer_parent(
    expr = fs::dir_delete(testdir),
    priority = "first"
  )

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
    pldsc = TRUE,
    name = "bipolar2021"
  )

  paths

}




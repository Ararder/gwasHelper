test_that("filepath_manager returns the correct filepath structure", {
  structure <- filepath_manager(test_path("data/gwas_repo/bipolar_disorder"), dir =test_path("data/gwas_repo/"))

  expected_structure <- list(
    base = fs::path(test_path("data/gwas_repo/bipolar_disorder")),
    base_job = fs::path(test_path("data/gwas_repo/bipolar_disorder/clean_bipolar_disorder.sh")),
    dataset_name = "bipolar_disorder",
    metafile = fs::path(test_path("data/gwas_repo/bipolar_disorder/meta.txt")),
    raw = fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder")),
    cleaned_dir = fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/")),
    clean = fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/cleaned_GRCh38.gz")),
    cleaned_meta = fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/cleaned_metadata.yaml")),
    sig_snps = fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/sig_snps.tsv")),
    analysis_dir = fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/analysis")),
    ldsc_dir = fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/analysis/ldsc")),
    ldsc_sumstats = fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/analysis/ldsc/ldsc")),
    ldsc_out  = fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/analysis/ldsc/ldsc_h2")),
    sbayes_dir =   fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/analysis", "sbayesr")),
    sbayes_ma =   fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/analysis", "sbayesr", "cleaned.ma")),
    sbayes_slurm_path =   fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/analysis", "sbayesr", "sbayes_bipolar_disorder.sh")),
    clumping_dir = fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/analysis", "clumping"))
  )

  expect_equal(structure, expected_structure)
})

# test_that("clean_sumstats_job writes a meta.txt file to the directory", {
#   tempdir <- withr::local_tempdir()
#   col_map <- c(
#     "col_SNP: SNP",
#     "col_POS: BP",
#     "col_BETA: Effect",
#     "col_SE: SE",
#     "col_P: Pval"
#   )
#   make_cleansumstats_job(tempdir, model = "logistic", col_map = col_map)
#
#   expect_true(
#     fs::file_exists(fs::path(tempdir, "meta.txt"))
#   )
# })

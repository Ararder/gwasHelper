test_that("filepath_manager returns the correct filepath structure", {
  structure <- filepath_manager(test_path("data/gwas_repo/bipolar_disorder"), dir =test_path("data/gwas_repo/"))

  expected_structure <- list(
    base = fs::path(test_path("data/gwas_repo/bipolar_disorder")),
    dataset_name = "bipolar_disorder",
    metafile = fs::path(test_path("data/gwas_repo/bipolar_disorder/meta.txt")),
    cleaned_dir = fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/")),
    clean = fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/cleaned_GRCh38.gz")),
    cleaned_meta = fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/cleaned_metadata.yaml")),
    sig_snps = fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/sig_snps.tsv")),
    analysis_dir = fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/analysis")),
    ldsc_dir = fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/analysis/ldsc")),
    ldsc_sumstats = fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/analysis/ldsc/ldsc")),
    ldsc_out  = fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/analysis/ldsc/ldsc_h2")),
    sbayes_dir =   fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/analysis", "sbayesr")),
    clumping_dir = fs::path(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/analysis", "clumping"))
  )

  expect_equal(structure, expected_structure)
})

# fs::file_touch(test_path("data/gwas_repo/bipolar_disorder/bipolar_disorder/cleaned_GRCh38.gz"))

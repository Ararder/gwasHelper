
### Parsers for LDSC
# ------------------------------------------------------------------------------


## first for the h2.log file
na_df <- dplyr::tibble(dataset_name=NA_character_, obs_h2=NA_real_,
              obs_se=NA_real_, lambda=NA_real_,
              mean_chi2=NA_real_, intercept=NA_real_,
              intercept_se=NA_real_, ratio=NA_real_)

test_that("Correctly parses ldsc file", {
  correct_df <- dplyr::tibble(dataset_name = "ldsc", obs_h2  = 0.4247, obs_se = 0.0446, lambda = 1.2564,
                mean_chi2 = 1.4588, intercept = 1.1155, intercept_se = 0.0099,
                ratio = 0.2518)

  expect_equal(parse_ldsc_h2(test_path("data/ldsc/h2/ldsc_h2.log")), correct_df)

})

test_that("If LDSC output file looks corrupt, then fails", {
  expect_equal(parse_ldsc_h2(test_path("data/ldsc/h2/ldsc_h2_err.log")), na_df)
})

# then output of --rg
test_that("If LDSC rg output is non-standard, return NA tibble",  {
  na_df <- dplyr::tibble(
    pheno1 = NA_character_, pheno2 = NA_character_, rg=NA_real_, rg_se=NA_real_, p = NA_real_)

  expect_equal(parse_ldsc_rg(test_path("data/ldsc/rg/antidep_perc_improv__heart_failure.log")), na_df)

})

parse_ldsc_rg(test_path("data/ldsc/rg/lymph__scz2022_eur.log"))


test_that("Correctly reads in LDSC rg results",  {
  real_df <- dplyr::tibble(
    pheno1 = "lymph", pheno2 = "scz2022_eur", rg=0.0807, se=0.0323,
    z = 2.4975, p = 0.0125, h2_obs = 0.3736, h2_obs_se = 0.016, h2_int = 1.0696,
    h2_int_se = 0.0163, gcov_int = -0.0274, gcov_int_se = 0.0156)


  expect_equal(parse_ldsc_rg(test_path("data/ldsc/rg/lymph__scz2022_eur.log")), real_df)

})


###  Parsers for SbayesR commandline
# ------------------------------------------------------------------------------

test_that("Parsing sbayes logs correctly", {
  real_df <- dplyr::tibble(
    dataset_name = "data", sbayes_h2 = 0.222098, sbayes_h2_se=0.026923)

  expect_equal(parse_sbayes_parres(test_path("data/sbayesr/adhd2019.parRes")), real_df)
})


test_that("SbayesR fails on non-standard file", {
  fail_df <- dplyr::tibble(
    dataset_name = NA_character_, sbayes_h2 = NA_real_, sbayes_h2_se = NA_real_)

  expect_equal(parse_sbayes_parres(test_path("data/sbayesr/corrupted.parRes")), fail_df)
})


df <- dplyr::tibble(data.table::fread("/nas/depts/007/sullilab/shared/gwas_sumstats/run2/mdd2019/mdd2019/analysis/magma/superclusters.gene_set.gsa.out", skip=3))
test_that("magma geneset parses correctly", {

  expect_true(length(colnames(parse_magma_geneset(test_path("data", "magma_geneset.tsv")))) == 7)

})




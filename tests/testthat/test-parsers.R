na_df <- dplyr::tibble(dataset_name=NA_character_, obs_h2=NA_character_,
              obs_se=NA_character_, lambda=NA_character_,
              mean_chi2=NA_character_, intercept=NA_character_,
              intercept_se=NA_character_, ratio=NA_character_)

test_that("Correctly parses ldsc file", {
  correct_df <- dplyr::tibble(dataset_name = ".", obs_h2  = 0.4247, obs_se = 0.0446, lambda = 1.2564,
                mean_chi2 = 1.4588, intercept = 1.1155, intercept_se = 0.0099,
                ratio = 0.2518)

  expect_equal(parse_ldsc_h2(test_path("data/ldsc_h2.log")), correct_df)

})

test_that("If LDSC output file looks corrupt, then fails", {
  expect_equal(parse_ldsc_h2(test_path("data/ldsc_h2_err.log")), na_df)
})

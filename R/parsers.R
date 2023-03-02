utils::globalVariables(c(".", "category"))

parse_ldsc_h2 <- function(path) {

  dataset_name <- get_analysis_phenotype(path)
  df <- readLines(path)

  if(length(df) != 32){
    return(dplyr::tibble(dataset_name=NA_character_, obs_h2=NA_real_,
                         obs_se=NA_real_, lambda=NA_real_,
                         mean_chi2=NA_real_, intercept=NA_real_,
                         intercept_se=NA_real_, ratio=NA_real_))
  }

  obs_h2 <- as.numeric(stringr::str_extract(df[26], "\\d{1}\\.\\d{1,5}"))

  obs_se <- stringr::str_extract(df[26], "\\(\\d{1}\\.\\d{1,5}") %>%
    stringr::str_remove(., "\\(") %>%
    as.numeric()

  lambda <- stringr::str_extract(df[27], " \\d{1}\\.\\d{1,5}") %>%
    as.numeric()

  mean_chi2 <- stringr::str_extract(df[28], " \\d{1}\\.\\d{1,5}") %>%
    as.numeric()

  intercept <- stringr::str_extract(df[29], " \\d{1}\\.\\d{1,5}") %>%
    as.numeric()

  intercept_se <- stringr::str_extract(df[29], "\\(\\d{1}\\.\\d{1,5}") %>%
    stringr::str_remove(., "\\(") %>%
    as.numeric()
  ratio <- stringr::str_extract(df[30], " \\d{1}\\.\\d{1,5}") %>%
    as.numeric()

  dplyr::tibble(dataset_name, obs_h2, obs_se, lambda, mean_chi2, intercept, intercept_se, ratio)
}

#' Title
#'
#' @param path path to LDSC rg output file
#'
#' @return a tibble
#' @export
#'
#' @examples \dontrun{
#' parse_ldsc_rg("path/to/ldsc_rg.results")
#' }
#'
parse_ldsc_rg <- function(path){
  strings <- readLines(path)

  if(length(strings) != 65) {
    return(dplyr::tibble(pheno1 = NA_character_, pheno2 = NA_character_, rg=NA_real_, rg_se=NA_real_, p = NA_real_))
  }


  pheno1 <- fs::path_file(path) %>%
    stringr::str_remove(".log") %>%
    stringr::str_split("__") %>% .[[1]] %>%
    .[1]

  pheno2 <- fs::path_file(path) %>%
    stringr::str_remove(".log") %>%
    stringr::str_split("__") %>% .[[1]] %>%
    .[2]

  names <- strings[61] %>%
    stringr::str_split(" ") %>%
    .[[1]] %>%
    .[. != ""] %>%
    .[-c(1:2)]

  vals <- strings[62] %>%
    stringr::str_split(" ") %>%
    .[[1]] %>%
    .[. != ""] %>%
    .[-c(1:2)] %>%
    as.numeric()

  dplyr::tibble(names, vals) %>%
    tidyr::pivot_wider(names_from = names, values_from = vals)  %>%
    dplyr::mutate(
      pheno1 = pheno1,
      pheno2 = pheno2
    ) %>%
    dplyr::select(pheno1, pheno2, dplyr::everything())

}

parse_sbayes_parres <- function(path){
  file <- suppressMessages(readr::read_tsv(path))

  if(length(readLines(path)) != 11) {
    message("non-standard file format. Returning NA")
    return(dplyr::tibble(dataset_name = NA_character_, sbayes_h2 = NA_real_, sbayes_h2_se = NA_real_))
  }

  res <- file %>%
    dplyr::slice(9) %>%
    dplyr::pull(1) %>%
    stringr::str_extract_all("\\d{1}\\.\\d{1,6}") %>%
    .[[1]]
  dplyr::tibble(dataset_name = get_analysis_phenotype(path), sbayes_h2 = as.numeric(res[1]), sbayes_h2_se = as.numeric(res[2]))
}

parse_clumping <- function(path) {
  dplyr::tibble(
    dataset_name = get_analysis_phenotype(path),
    n_loci = length(readLines(path))
  )
}

parse_sig_snps <- function(path){
  dataset_name = get_analysis_phenotype(path)
  suppressMessages(readr::read_tsv(path)) %>%
    dplyr::mutate(dataset_name = dataset_name) %>%
    dplyr::select(2,1)
}

parse_pldsc <- function(path) {
  name <- janitor::make_clean_names(
    fs::path_file(path) %>% stringr::str_remove(".results")
  )
  suppressMessages(readr::read_tsv(path)) |>
    janitor::clean_names() |>
    dplyr::filter(category == "L2_1") |>
    dplyr::mutate(category =  name)
}

parse_magma_geneset <- function(path){
  dplyr::tibble(data.table::fread(path, skip = 3)) |>
    dplyr::mutate(VARIABLE = janitor::make_clean_names(FULL_NAME)) |>
    dplyr::select(-FULL_NAME)
}

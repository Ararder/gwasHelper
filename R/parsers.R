utils::globalVariables(c("."))
# parses the output of a standard ldsc commandline call to ldsc.py --h2
# using regex pattern matching
  # extract the output of LDSC h2
parse_ldsc_h2 <- function(path) {

  dataset_name <- get_analysis_phenotype(path)
  df <- readLines(path)

  if(length(df) != 32){
    return(dplyr::tibble(dataset_name=NA_character_, obs_h2=NA_character_,
                         obs_se=NA_character_, lambda=NA_character_,
                         mean_chi2=NA_character_, intercept=NA_character_,
                         intercept_se=NA_character_, ratio=NA_character_))
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

# get_sbayes_h2 <- function(path){
#   res <- read_tsv(path) %>%
#     slice(9) %>%
#     pull(1) %>%
#     str_extract_all("\\d{1}\\.\\d{1,6}") %>%
#     .[[1]]
#   tibble(dataset_name = get_analysis_phenotype(path), sbayes_h2 = res[1], sbayes_h2_se = res[2])
# }
#
#
# collect_sbayes <- function(dir = "/nas/depts/007/sullilab/shared/gwas_sumstats/run2"){
#   map_df(dir_ls(dir, glob ="*parRes", recurse=TRUE), get_sbayes_h2) %>%
#     mutate(across(c(2,3), as.double))
# }
#
# collect_ldsc <- function(dir = "/nas/depts/007/sullilab/shared/gwas_sumstats/run2"){
#   map_df(dir_ls(dir, glob ="*ldsc_h2.log", recurse=TRUE), get_ldsc_h2)
# }
#
#
# get_cleaned_sumstats <- function(dir = "/nas/depts/007/sullilab/shared/gwas_sumstats/run2"){
#   tibble(clean = dir_ls(dir, glob ="*cleaned_GRCh38.gz", recurse=TRUE)) %>%
#     dplyr::mutate(dataset_name = path_file(path_dir(clean))) %>%
#     dplyr::select(2,1)
# }
#
# dir_ls("run2", glob = "*sig_snps.tsv", recurse=TRUE)
# get_sig_snps <- function(path){
#   dataset_name = get_analysis_phenotype(path)
#   read_tsv(path) %>%
#     mutate(dataset_name = dataset_name) %>%
#     select(2,1)
# }
# collect_sig_snps <- function(dir = "/nas/depts/007/sullilab/shared/gwas_sumstats/run2"){
#   map_df(dir_ls(dir, glob = "*sig_snps.tsv", recurse=TRUE), get_sig_snps)
# }
#
#
#
# collect_loci <- function(dir = "/nas/depts/007/sullilab/shared/gwas_sumstats/run2"){
#   get_loci <- function(path) {
#     tibble(
#       dataset_name = get_analysis_phenotype(path),
#       n_loci = length(readLines(path))
#     )
#
#   }
#   map_df(dir_ls(dir, glob = "*n_loci.bedtools", recurse=TRUE), get_loci)
# }
#
#
#
# unique_order <- function(string1, string2) {
#   string <- sort(c(string1, string2))
#   stringr::str_c(string[1], string[2], sep = "__")
#
# }
# unique_order <- Vectorize(unique_order)

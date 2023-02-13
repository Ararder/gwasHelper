utils::globalVariables(c("sumstat_path"))


collect_sbayes <- function(dir = Sys.getenv("gwasHelper_repo")){
  purrr::map_df(fs::dir_ls(dir, glob ="*parRes", recurse=TRUE), parse_sbayes_parres) %>%
    dplyr::mutate(dplyr::across(c(2,3), as.double))
}


collect_ldsc <- function(dir = Sys.getenv("gwasHelper_repo")){
  purrr::map_df(fs::dir_ls(dir, glob ="*ldsc_h2.log", recurse=TRUE), parse_ldsc_h2)
}

collect_sumstats <- function(dir = Sys.getenv("gwasHelper_repo")){

    dplyr::tibble(sumstat_path = fs::dir_ls(dir, glob ="*cleaned_GRCh38.gz", recurse=TRUE)) %>%
    dplyr::mutate(dataset_name = fs::path_file(fs::path_dir(sumstat_path))) %>%
    dplyr::select(2,1)
}


collect_sig_snps <- function(dir = Sys.getenv("gwasHelper_repo")){
  purrr::map_df(fs::dir_ls(dir, glob = "*sig_snps.tsv", recurse=TRUE), parse_sig_snps)
}

collect_loci <- function(dir = Sys.getenv("gwasHelper_repo")){
  purrr::map_df(fs::dir_ls(dir, glob = "*n_loci.bedtools", recurse=TRUE), parse_clumping)
}



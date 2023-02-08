utils::globalVariables(c("clean"))
#' Get filepath and dataset_name for all cleaned sumstats
#'
#' @param dir filepath for the directory to search. Defaults to use the ENV
#' variable "gwasHelper_repo"
#'
#' @return a tibble
#' @export
#'
#' @examples \dontrun{
#' # most common use is to call without arugemnt
#' get_cleaned_sumstats()
#' }
get_cleaned_sumstats <- function(dir = Sys.getenv("gwasHelper_repo")){
  dplyr::tibble(clean = fs::dir_ls(dir, glob ="*cleaned_GRCh38.gz", recurse=TRUE)) %>%
    dplyr::mutate(dataset_name = fs::path_file(fs::path_dir(clean))) %>%
    dplyr::select(2,1)
}

get_analysis_phenotype <- function(path){
  # get the phenotype from a file in analysis
  fs::path_dir(path) %>%
    fs::path_dir() %>%
    stringr::str_remove("/analysis") %>%
    fs::path_file()
}



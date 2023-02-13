get_analysis_phenotype <- function(path){
  # get the phenotype from a file in analysis
  fs::path_dir(path) %>%
    fs::path_dir() %>%
    stringr::str_remove("/analysis") %>%
    fs::path_file()
}

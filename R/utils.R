#' Title
#'
#' @param time --time parameter in slurm
#' @param mem  --mem parameter in slurm
#' @param output_dir in what directory to put the slurm output file
#' @param dependency default null, if specified will add dependecy=afterok:{dependency}
#'
#' @return a character vector
#' @export
#'
#' @examples slurm_header(24, 13, "my_dir/run_job")
slurm_header <- function(time,mem, output_dir, dependency=NULL) {
  c(
    "#!/bin/bash",
    glue::glue("#SBATCH --time={time}:00:00"),
    glue::glue("#SBATCH --mem={mem}gb"),
    glue::glue("#SBATCH --output={output_dir}/slurm-%j.out"),
    glue::glue("#SBATCH --dependency=afterok:{dependency}"),
    "\n"
  )
}
slurm_header <- Vectorize(slurm_header)

get_analysis_phenotype <- function(path){
  # get the phenotype from a file in analysis
  fs::path_dir(path) %>%
    fs::path_dir() %>%
    stringr::str_remove("/analysis") %>%
    fs::path_file()
}

validate_columns <- function(type, df) {
  if(type == "magma") {
      stopifnot("Magma requires RSID, POS, CHR, N and P columns" =all(c("RSID", "POS", "CHR", "P", "N") %in% colnames(df)))
    }

}

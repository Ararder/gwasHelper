slurm_header <- function(time,mem, output, dependency=NULL) {
  c(
    "#!/bin/bash",
    glue::glue("#SBATCH --time={time}:00:00"),
    glue::glue("#SBATCH --mem={mem}gb"),
    glue::glue("#SBATCH --output={output}/slurm-%j.out"),
    glue::glue("#SBATCH --dependency=afterok:{dependency}"),
    "\n"
  ) %>%
    list() %>%
    purrr::reduce(c)
}
slurm_header <- Vectorize(slurm_header)

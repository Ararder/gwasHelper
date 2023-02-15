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

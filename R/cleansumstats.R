clean_sumstats <- function(infile, name,model = "logistic", col_map) {
  if(missing(name)) {
    paths <- filepath_manager(infile)
  } else {
    paths <- filepath_manager(infile, name)
  }


  # create the new folder, and copy over the file
  fs::file_copy(infile, fs::dir_create(paths[["base"]]))

  # create the slurm header for the job
  slurm_header <- create_slurm_header(paths[["base"]])

  # create the cleansumstats job --> will also writeout metadata file
  cleansumstats_job <- make_cleansumstats_job(
    paths[["base"]],
    model = model,
    col_map = col_map
  )

  # initialize output folders for secondary analysis
  create_dirs <- glue::glue(
    "mkdir -p {paths[['ldsc_dir']]}",
    " {paths[['clumping_dir']]}",

  )
  # LDSC
  ldsc_job <- run_ldsc(paths)

  # Clumping

  # SbayesR (?)

  c(slurm_header, cleansumstats_job, create_dirs, ldsc_job)
}

make_cleansumstats_job <- function(dir, model, col_map) {
  # setup and create the metadata file, and write it to the directory
  header <- c(
    "cleansumstats_metafile_kind: minimal",
    glue::glue("path_sumStats: {fs::path_file(dir)}"),
    glue::glue("stats_Model: {model}")
  )
  # col_map is passed
  writeLines(c(header, col_map), fs::path(dir, "meta.txt"))

  glue::glue(
    "{Sys.getenv('cleanSumstats')}/cleansumstats.sh ",
    "-i {fs::path(dir, 'meta.txt')} ",
    "-d {Sys.getenv('cleanSumstats')}/out_dbsnp ",
    "-k {Sys.getenv('cleanSumstats')}/out_1kgp ",
    "-o {dir}/{fs::path_file(dir)}",
  )
}
create_slurm_header <- function(dir) {
  c(
    "#!/bin/bash",
    "#SBATCH --time=3:00:00",
    "#SBATCH --mem=10g",
    glue::glue("#SBATCH --output={dir}/slurm-%j.out")
  )
}

clean_sumstats <- function(infile, col_map, name,model = "logistic") {
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
  clump_job <- run_clumping(paths[["clean"]], paths[["clumping_dir"]])

  # calculate sig sns
  sig_snps <- glue::glue(
  "R -e ",
  "'readr::write_tsv(dplyr::tibble(n_sig = nrow(dplyr::filter(data.table::fread(commandArgs(trailingOnly=TRUE)[1]), P < 5e-08))),commandArgs(trailingOnly=TRUE)[2])'",
  " --args {paths[['clean']]} {paths[['sig_snps']]}"
  )

  cleanup <- glue::glue(
    "rm {paths[['raw']]}"
  )
  # SbayesR (?)

  c(slurm_header, cleansumstats_job, create_dirs, ldsc_job,
    clump_job, sig_snps, cleanup)
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

filepath_manager <- function(infile, dir=Sys.getenv("gwasHelper_repo"), name) {
  if(missing(name)) {
    name <- fs::path_ext_remove(fs::path_ext_remove(fs::path_ext_remove(fs::path_file(infile))))
  }
  base <- fs::path(dir, name)

  # Base level
  dataset_name <- fs::path_file(base)
  metafile <-     fs::path(base, "meta.txt")
  raw <- fs::path(base, fs::path_file(infile))
  base_job <- fs::path(base, paste0("clean", "_", dataset_name, ".sh"))

  # dir storing cleaned GWAS + post analysis
  cleaned_dir <-  fs::path(base, dataset_name)
  clean =         fs::path(cleaned_dir, "cleaned_GRCh38.gz")
  cleaned_meta <- fs::path(cleaned_dir, "cleaned_metadata.yaml")
  sig_snps <-     fs::path(cleaned_dir, "sig_snps.tsv")

  # dir with different post GWAS analysis
  analysis_dir <- fs::path(cleaned_dir, "analysis")

  ldsc_dir <-     fs::path(analysis_dir, "ldsc")
  ldsc_sumstats <-fs::path(analysis_dir, "ldsc", "ldsc")
  ldsc_out <-fs::path(analysis_dir, "ldsc", "ldsc_h2")

  sbayes_dir <-   fs::path(analysis_dir, "sbayesr")
  sbayes_ma <-   fs::path(analysis_dir, "sbayesr", "cleaned.ma")
  sbayes_slurm_path = fs::path(analysis_dir, "sbayesr", paste0("sbayes", "_", paths[["dataset_name"]], ".sh"))
  clumping_dir <- fs::path(analysis_dir, "clumping")


  list(
    "base" = base,
    "base_job" =  base_job,
    "dataset_name" = dataset_name,
    "metafile" = metafile,
    "raw" = raw,

    "cleaned_dir" = cleaned_dir,
    "clean" = clean,
    "cleaned_meta" = cleaned_meta,
    "sig_snps" = sig_snps,

    "analysis_dir" = analysis_dir,
    "ldsc_dir" = ldsc_dir,
    "ldsc_sumstats" = ldsc_sumstats,
    "ldsc_out" = ldsc_out,
    "sbayes_dir" = sbayes_dir,
    "sbayes_ma" = sbayes_ma,
    "sbayes_slurm_path" = sbayes_slurm_path,
    "clumping_dir" = clumping_dir

  )

}


match_cols <- function(
  col_CHR = NULL, col_POS = NULL, col_BETA = NULL,
  col_SNP = NULL, col_EffectAllele = NULL,
  col_SE = NULL, col_Z = NULL, col_OR = NULL,
  col_N = NULL, col_CaseN = NULL, col_ControlN = NULL,
  col_EAF = NULL, col_INFO = NULL, col_StudyN=NULL,
  col_OtherAllele = NULL, col_P = NULL,
  stats_CaseN = NULL, stats_ControlN=NULL,
  stats_StudyN=NULL) {

  list(
    glue::glue("col_CHR: {col_CHR}"),
    glue::glue("col_POS: {col_POS}"),
    glue::glue("col_SNP: {col_SNP}"),
    glue::glue("col_BETA: {col_BETA}"),
    glue::glue("col_EffectAllele: {col_EffectAllele}"),
    glue::glue("col_SE: {col_SE}"),
    glue::glue("col_Z: {col_Z}"),
    glue::glue("col_OR: {col_OR}"),
    glue::glue("col_N: {col_N}"),
    glue::glue("col_P: {col_P}"),
    glue::glue("col_CaseN: {col_CaseN}"),
    glue::glue("col_ControlN: {col_ControlN}"),
    glue::glue("col_EAF: {col_EAF}"),
    glue::glue("col_INFO: {col_INFO}"),
    glue::glue("col_StudyN: {col_StudyN}"),
    glue::glue("col_OtherAllele: {col_OtherAllele}"),
    glue::glue("stats_CaseN: {stats_CaseN}"),
    glue::glue("stats_ControlN: {stats_ControlN}"),
    glue::glue("stats_StudyN: {stats_StudyN}")
  ) %>%
    purrr::reduce(c)
}

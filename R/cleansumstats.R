#' Title
#'
#' @param infile raw gwas sumstat to clean
#' @param col_map map of column names
#' @param name gwas sumstat name
#' @param model logistic/linear
#'
#' @return a filepath to the slurm job
#' @export
#'
#' @examples \dontrun{
#' clean_sumstats("my_raw_gwas.tsv.gz",name =  "mdd2024",col_map = col_map)
#' }
clean_sumstats <- function(infile, col_map, name, model = "logistic") {
  if(missing(name)) {
    paths <- filepath_manager(infile)

    sbayes_job <- glue::glue("Rscript -e 'gwasHelper::sbayes_run(gwasHelper::filepath_manager(commandArgs(trailingOnly=TRUE)))'") %>%
      glue::glue(" --args {infile}")
  } else {
    paths <- filepath_manager(infile, name)

    sbayes_job <- glue::glue("Rscript -e 'gwasHelper::sbayes_run(gwasHelper::filepath_manager(commandArgs(trailingOnly=TRUE)[2],commandArgs(trailingOnly=TRUE)[3]))'") %>%
      glue::glue(" --args {infile} {name}")
  }



  # create the new folder, and copy over the file
  fs::file_copy(infile, fs::path(fs::dir_create(paths[["base"]]), fs::path_file(infile)))

  # create the slurm header for the job
  slurm_header <- slurm_header(time=3, mem = 4, output_dir = paths[["base"]])

  # create the cleansumstats job --> will also writeout metadata file
  cleansumstats_job <- make_cleansumstats_job(
    paths,
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
  "Rscript -e ",
  "'readr::write_tsv(dplyr::tibble(n_sig = nrow(dplyr::filter(data.table::fread(commandArgs(trailingOnly=TRUE)[2]), P < 5e-08))),commandArgs(trailingOnly=TRUE)[3])'",
  " --args {paths[['clean']]} {paths[['sig_snps']]}"
  )

  cleanup <- glue::glue(
    "rm {paths[['raw']]}"
  )
  # code to sbatch sbayes
  sbatch_sbayes <- glue::glue("sbatch {paths[['sbayes_slurm_path']]}")

  captured <- c(slurm_header, cleansumstats_job, create_dirs, ldsc_job,
    clump_job,sbayes_job,sbatch_sbayes, sig_snps, cleanup)

  # run s-ldsc?
  fs::dir_create(paths[["pldsc_siletti"]], recurse=TRUE)
  writeLines(
    c(
      cleansumstats_pldsc(paths, Sys.getenv("siletti_clusters_ldscores")),
      "\n",
      cleansumstats_pldsc(paths, Sys.getenv("siletti_superclusters_ldscores"))
    ),
    paths[['pldsc_siletti_slurm']]
  )
  sldsc_siletti <- glue::glue("chmod 700 {paths[['pldsc_siletti_slurm']]} && {paths[['pldsc_siletti_slurm']]}")


  captured <- c(
    slurm_header,
    cleansumstats_job,
    create_dirs,
    ldsc_job,
    clump_job,
    sldsc_siletti,
    sbayes_job,
    sbatch_sbayes,
    sig_snps,
    cleanup)

  writeLines(captured, paths[["base_job"]])

  paths[["base_job"]]
}


make_cleansumstats_job <- function(paths, model, col_map) {
  # setup and create the metadata file, and write it to the directory
  header <- c(
    "cleansumstats_metafile_kind: minimal",
    glue::glue("path_sumStats: {fs::path_file(paths[['raw']])}"),
    glue::glue("stats_Model: {model}")
  )
  # col_map is passed
  writeLines(c(header, col_map), paths[['metafile']])
  paths[['metafile']]
  glue::glue(
    "{Sys.getenv('cleanSumstats')}/cleansumstats.sh ",
    "-i {paths[['metafile']]} ",
    "-d {Sys.getenv('cleanSumstats')}/out_dbsnp ",
    "-k {Sys.getenv('cleanSumstats')}/out_1kgp ",
    "-o {paths[['cleaned_dir']]}",
  )
}



#' Title
#'
#' @param infile gwas sumstat to clean
#' @param dir directory used to store cleaned sumstats
#' @param name optional: name to call the sumstat
#'
#' @return a list of filepaths
#' @export
#'
#' @examples \dontrun{
#' path <- filepath_manager("adhd_2023.tsv.gz", "adhd_pgc2023")
#' }
filepath_manager <- function(infile, name,dir) {
  if(missing(dir)) {
    dir <- Sys.getenv("gwasHelper_repo")
  }
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
  pldsc <- fs::path(ldsc_dir, "pldsc")
  pldsc_siletti <- fs::path(pldsc, "siletti")
  pldsc_siletti_slurm <- fs::path(pldsc, "siletti", "run_all.sh")

  sbayes_dir <-   fs::path(analysis_dir, "sbayesr")
  sbayes_ma <-   fs::path(analysis_dir, "sbayesr", "cleaned.ma")
  sbayes_slurm_path = fs::path(analysis_dir, "sbayesr", paste0("sbayes", "_", dataset_name, ".sh"))
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
    "pldsc" = pldsc,
    "pldsc_siletti" = pldsc_siletti,
    "pldsc_siletti_slurm" = pldsc_siletti_slurm,
    "sbayes_dir" = sbayes_dir,
    "sbayes_ma" = sbayes_ma,
    "sbayes_slurm_path" = sbayes_slurm_path,
    "clumping_dir" = clumping_dir

  )

}


#' Title
#'
#' @param col_CHR chromosome
#' @param col_POS position
#' @param col_BETA effect column (beta or odds ratio)
#' @param col_SNP rsid column
#' @param col_EffectAllele effect allele
#' @param col_SE standard error
#' @param col_Z Z score
#' @param col_OR odds ratio
#' @param col_N sample size
#' @param col_CaseN case sample size
#' @param col_ControlN control sample size
#' @param col_EAF allele frequency of effect allele
#' @param col_INFO info column
#' @param col_StudyN development column, dont use
#' @param col_OtherAllele other allele
#' @param col_P P value
#' @param stats_CaseN case sample size for entire study(not per snp)
#' @param stats_ControlN control sample size for entire study (not per snp)
#' @param stats_StudyN (sample size for entire study - not per snp)
#'
#' @return a col_map to use with clean_sumstats
#' @export
#'
#' @examples \dontrun{
#' col_map(col_CHR = "CHR")
#' }
#'
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

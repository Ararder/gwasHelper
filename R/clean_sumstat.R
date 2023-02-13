

filepath_manager <- function(infile, dir=Sys.getenv("gwasHelper_repo"), name) {
  if(missing(name)) {
    name <- fs::path_ext_remove(fs::path_ext_remove(fs::path_ext_remove(fs::path_file(infile))))
  }
  base <- fs::path(dir, name)

  # Base level
  dataset_name <- fs::path_file(base)
  metafile <-     fs::path(base, "meta.txt")

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
  clumping_dir <- fs::path(analysis_dir, "clumping")


  list(
    "base" = base,
    "dataset_name" = dataset_name,
    "metafile" = metafile,

    "cleaned_dir" = cleaned_dir,
    "clean" = clean,
    "cleaned_meta" = cleaned_meta,
    "sig_snps" = sig_snps,

    "analysis_dir" = analysis_dir,
    "ldsc_dir" = ldsc_dir,
    "ldsc_sumstats" = ldsc_sumstats,
    "ldsc_out" = ldsc_out,
    "sbayes_dir" = sbayes_dir,
    "clumping_dir" = clumping_dir

  )

}

#
#
#
#
# run_gwas_cleaning <- function(path, col_map, model="logistic",
#                               name,
#                               test_folder="/nas/depts/007/sullilab/shared/gwas_sumstats/run2",
#                               sbayes=TRUE) {
#
#   # 1. Setup output direcgtory
#   ## ----------------------------------------------------------------------------
#   ## remove 3x extensions from filepath,
#   ## copy over gwas file and metadata to all be in the same directory
#   ## create the wrapper for calling cleansumstats.sh
#
#   if(missing(name)) name <- path_ext_remove(path_ext_remove(path_ext_remove(fs::path_file(path))))
#   dir <- dir_create(here(path(test_folder, name)))
#   analysis_dir <- path(dir, name, "analysis")
#   ldsc_dir <- path(analysis_dir, "ldsc")
#   file_copy(path, dir)
#   sbatch_out <- paste0(dir, "/run.sh")
#   sbatch_job <- make_cleansumstats_job(dir)
#
#   # define the slurm header
#   slurm_header <- c(
#     "#!/bin/bash",
#     "#SBATCH --time=3:00:00",
#     "#SBATCH --mem=10g",
#     glue("#SBATCH --output={dir}/slurm-%j.out")
#   )
#
#
#   # 2. create the metadata file
#   ##  ----------------------------------------------------------------------------
#   header <- c(
#     "cleansumstats_metafile_kind: minimal",
#     glue("path_sumStats: {fs::path_file(path)}"),
#     glue("stats_Model: {model}")
#   )
#
#   # create filepath for metadata, and write it out
#   writeLines(c(header, col_map), paste0(dir, "/", "meta.txt"))
#
#   # 3. Create LDSC call add  to slurm job
#   ##  ------------------------------------------------------------------------------
#   sumstat_path <- path(dir, name, "cleaned_GRCh38.gz")
#
#   ldsc <- c("module unload python/3.9.6", "module load ldsc", run_ldsc(sumstat_path,col_map))
#
#   # Cleanup
#   raw_sumstats <- path(dir, name)
#   remove_raw <- glue::glue("rm {path(dir, fs::path_file(path))}")
#   analysis_dir <- glue::glue("mkdir -p {analysis_dir}/ldsc")
#   cleanup <- c(remove_raw, analysis_dir)
#
#   # 4. Add a Rscript to check n signifcant snps
#   # ------------------------------------------------------------------------------
#
#   job4 <- glue("Rscript /nas/depts/007/sullilab/shared/gwas_sumstats/R/scripts/eval_sig_snps.R {sumstat_path}")
#
#   # 5. Add all the different parts together, then sbatch job (check sbayes)
#   ## ------------------------------------------------------------------------------
#   if(sbayes){
#     sbayes_job <- glue("Rscript /nas/depts/007/sullilab/shared/gwas_sumstats/R/scripts/sbayes_cleaning_pipeline.R {sumstat_path}")
#     merged <- c(slurm_header, "\n", sbatch_job, "\n",cleanup,"\n", sbayes_job, job4, "\n" ,ldsc, "\n")
#   } else {
#     merged <- c(slurm_header, "\n", sbatch_job, "\n",cleanup,"\n", job4, "\n" ,ldsc, "\n")
#
#   }
#
#
#   writeLines(merged, sbatch_out)
#   sbatch_out
#   # system(glue("sbatch {sbatch_out}"))
#
# }
#
# match_cols <- function(
#   col_CHR = NULL, col_POS = NULL, col_BETA = NULL,
#   col_SNP = NULL, col_EffectAllele = NULL,
#   col_SE = NULL, col_Z = NULL, col_OR = NULL,
#   col_N = NULL, col_CaseN = NULL, col_ControlN = NULL,
#   col_EAF = NULL, col_INFO = NULL, col_StudyN=NULL,
#   col_OtherAllele = NULL, col_P = NULL,
#   stats_CaseN = NULL, stats_ControlN=NULL,
#   stats_StudyN=NULL) {
#
#   list(
#     glue::glue("col_CHR: {col_CHR}"),
#     glue::glue("col_POS: {col_POS}"),
#     glue::glue("col_SNP: {col_SNP}"),
#     glue::glue("col_BETA: {col_BETA}"),
#     glue::glue("col_EffectAllele: {col_EffectAllele}"),
#     glue::glue("col_SE: {col_SE}"),
#     glue::glue("col_Z: {col_Z}"),
#     glue::glue("col_OR: {col_OR}"),
#     glue::glue("col_N: {col_N}"),
#     glue::glue("col_P: {col_P}"),
#     glue::glue("col_CaseN: {col_CaseN}"),
#     glue::glue("col_ControlN: {col_ControlN}"),
#     glue::glue("col_EAF: {col_EAF}"),
#     glue::glue("col_INFO: {col_INFO}"),
#     glue::glue("col_StudyN: {col_StudyN}"),
#     glue::glue("col_OtherAllele: {col_OtherAllele}"),
#     glue::glue("stats_CaseN: {stats_CaseN}"),
#     glue::glue("stats_ControlN: {stats_ControlN}"),
#     glue::glue("stats_StudyN: {stats_StudyN}")
#   ) %>%
#     reduce(c)
# }

#' Title
#'
#' @param infile raw gwas sumstat to clean
#' @param col_map map of column names
#' @param name gwas sumstat name
#' @param model logistic/linear
#' @param magma_gene run magma gene analysis?
#' @param magma_gene_set run magma gene analysis?
#' @param pldsc run partitioned ldscore regression?
#' @param sbayesr run sbayesR to generate a PGS?
#'
#'
#' @return a filepath to the slurm job
#' @export
#'
#' @examples \dontrun{
#' clean_sumstats("my_raw_gwas.tsv.gz",name =  "mdd2024",col_map = col_map)
#' }
clean_sumstats <- function(infile,
                           col_map,
                           name,
                           model = "logistic",
                           magma_gene = FALSE,
                           magma_gene_set = FALSE,
                           pldsc = FALSE,
                           sbayesr = FALSE) {
  if(missing(name)) {
    paths <- filepath_manager(infile)

    sbayes_job <- glue::glue("Rscript -e 'gwasHelper::sbayes_run(gwasHelper::filepath_manager(commandArgs(trailingOnly=TRUE)))'") %>%
      glue::glue(" --args {infile}")
  } else {
    paths <- filepath_manager(infile, name)

    sbayes_job <- glue::glue("Rscript -e 'gwasHelper::sbayes_run(gwasHelper::filepath_manager(commandArgs(trailingOnly=TRUE)[2],commandArgs(trailingOnly=TRUE)[3]))'") %>%
      glue::glue(" --args {infile} {name}")
  }


  ## Basic analysis of in the clean_sumstats pipeline
  # ----------------------------------------------------------------------------

  # create the new folder, and copy over the file
  fs::file_copy(infile, fs::path(fs::dir_create(paths[["base"]]), fs::path_file(infile)))

  # create the slurm header for the job
  slurm_header <- slurm_header(time=24, mem = 4, output_dir = paths[["base"]])

  # create the cleansumstats job --> will also writeout metadata file
  cleansumstats_job <- make_cleansumstats_job(
    paths,
    model = model,
    col_map = col_map
  )
  # Limit pvalues to minimum precision
  adjust <-correct_small_pval(paths)

  # initialize output folders for secondary analysis
  create_dirs <- initialize_folder_structure(paths)

  # LDSC, Clumping, number of significant SNPs and cleanup
  ldsc_job <- run_ldsc(paths)
  clump_job <- run_clumping(paths[["clean"]], paths[["clumping_dir"]])
  sig_snps <- n_significant_snps(paths)
  cleanup <- cleanup_files(paths)

  # Collect the code
  captured <- c(
    slurm_header,
    "# 1) Cleansumstats pipeline, to QC and harmonize summary statistic",
    cleansumstats_job,
    "\n",
    "# 2) If Pvalues are extremely small, set them to e-308",
    adjust,
    "# Create directories, and calculate the number of significant SNPs",
    create_dirs,
    sig_snps,
    "\n",
    "# 2) LDSC to munge and estimate heritability",
    ldsc_job,
    "\n",
    "# 3) Clumping to get number of loci",
    clump_job,
    "\n"
  )


  # P-LDSC
  if(pldsc) {

    ldscores <- c(
      Sys.getenv("siletti_superclusters_ldscores"),
      Sys.getenv("siletti_clusters_ldscores")
      )

    pldsc <- cleansumstats_pldsc(paths,ldscores)
    captured <-  c(
      captured,
      "# 4) Stratified LDSC using scRNA human brain cell data",
      pldsc,
      "\n "
    )
  }
  ## MAGMA
  # ----------------------------------------------------------------------------
  if(magma_gene) {
    magma <- ""
    captured <- c(
      captured,
      "# 5) Generate code to run Magma gene-set ",
      magma,
      "\n"
    )
  }

  ## SbayesR
  # ----------------------------------------------------------------------------

  if(sbayesr) {
    sbatch_sbayes <- glue::glue("sbatch {paths[['sbayes_slurm_path']]}")
    captured <- c(
      captured,
        "# 6) Generate code to run SbayesR for polygenic scores ",
        sbayes_job,
        sbatch_sbayes,
        "\n"
    )
  }

  ## cleanup code
  # ----------------------------------------------------------------------------
  captured <- c(
    captured,
    "# Clean up files that are no longer needed",
    cleanup
    )

  writeLines(captured, paths[["base_job"]])

  paths[["base_job"]]
}


initialize_folder_structure <- function(paths) {
  glue::glue("mkdir -p {paths[['ldsc_dir']]} {paths[['clumping_dir']]}")

}




n_significant_snps <- function(paths){
  glue::glue("gunzip -c {paths[['clean']]} | awk -F '\\t' 'NR==1 {{for (i=1; i<=NF; i++) {{if ($i==\"P\") {{col=i; break}}}}}} NR>1 && $col <= 5e-08' | wc -l > {paths[['sig_snps']]}")
}

cleanup_files <- function(paths) {
  # remove the raw sumstat
  glue::glue("rm {paths[['raw']]}")

}

make_cleansumstats_job <- function(paths, model, col_map) {
  # setup and create the metadata file, and write it to the directory
  header <- c(
    "cleansumstats_metafile_kind: minimal",
    glue::glue("path_sumStats: {fs::path_file(paths[['raw']])}"),
    glue::glue("stats_Model: {model}")
  )
  # Write out the metafile
  writeLines(c(header, col_map), paths[['metafile']])


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
  temp =         fs::path(cleaned_dir, "temp.gz")
  cleaned_meta <- fs::path(cleaned_dir, "cleaned_metadata.yaml")
  sig_snps <-     fs::path(cleaned_dir, "sig_snps.tsv")
  build37 <-      fs::path(cleaned_dir, "cleaned_GRCh37.gz")

  # dir with different post GWAS analysis
  analysis_dir <- fs::path(cleaned_dir, "analysis")

  # basic ldsc + partitioned ldsc
  ldsc_dir <-     fs::path(analysis_dir, "ldsc")
  ldsc_sumstats <-fs::path(analysis_dir, "ldsc", "ldsc")
  ldsc_out <-fs::path(analysis_dir, "ldsc", "ldsc_h2")
  pldsc <- fs::path(ldsc_dir, "pldsc")
  pldsc_siletti <- fs::path(pldsc, "siletti")
  pldsc_slurm <- fs::path(pldsc, "run_all.sh")

  # magma
  magma_dir <- fs::path(analysis_dir, "magma")

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
    "temp" = temp,
    "cleaned_meta" = cleaned_meta,
    "sig_snps" = sig_snps,
    "build37" = build37,

    "analysis_dir" = analysis_dir,
    "ldsc_dir" = ldsc_dir,
    "ldsc_sumstats" = ldsc_sumstats,
    "ldsc_out" = ldsc_out,
    "pldsc" = pldsc,
    "pldsc_siletti" = pldsc_siletti,
    "pldsc_slurm" = pldsc_slurm,

    "magma_dir"  = magma_dir,

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
#' match_cols(col_CHR = "CHR")
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
    glue::glue("col_N: '{col_N}'"),
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
#' Title
#'
#' @param infile raw gwas sumstat to clean
#' @param name gwas sumstat name
#' @param meta filepath to metafile
#'
#' @return a filepath to the slurm job
#' @export
#'
#' @examples \dontrun{
#' migrate("my_raw_gwas.tsv.gz",name =  "mdd2024",col_map = col_map)
#' }
migrate <- function(infile, meta, name) {
  if(missing(name)) {
    paths <- filepath_manager(infile)

    sbayes_job <- glue::glue("Rscript -e 'gwasHelper::sbayes_run(gwasHelper::filepath_manager(commandArgs(trailingOnly=TRUE)))'") %>%
      glue::glue(" --args {infile}")
  } else {
    # the name variable needs to be passed to SbayesR helper
    paths <- filepath_manager(infile, name)
    sbayes_job <- glue::glue("Rscript -e 'gwasHelper::sbayes_run(gwasHelper::filepath_manager(commandArgs(trailingOnly=TRUE)[2],commandArgs(trailingOnly=TRUE)[3]))'") %>%
      glue::glue(" --args {infile} {name}")
  }




  cleansumstats_job <- glue::glue(
    "{Sys.getenv('cleanSumstats')}/cleansumstats.sh ",
    "-i {paths[['metafile']]} ",
    "-d {Sys.getenv('cleanSumstats')}/out_dbsnp ",
    "-k {Sys.getenv('cleanSumstats')}/out_1kgp ",
    "-o {paths[['cleaned_dir']]} ",
    "-b /work/users/a/r/arvhar/tmp"
  )




  # create the new folder, and copy over the file
  fs::file_copy(infile, fs::path(fs::dir_create(paths[["base"]]), fs::path_file(infile)))

  # copy over metafile
  fs::file_copy(meta, paths[["metafile"]])


  # create the slurm header for the job
  slurm_header <- slurm_header(time=24, mem = 4, output_dir = paths[["base"]])


  # initialize output folders for secondary analysis
  create_dirs <- initialize_folder_structure(paths)

  # LDSC
  ldsc_job <- run_ldsc(paths)

  # Clumping
  clump_job <- run_clumping(paths[["clean"]], paths[["clumping_dir"]])

  # Find the nubmer of significant snps
  sig_snps <- n_significant_snps(paths)
  fix_pvals <-correct_small_pval(paths)

  cleanup <- cleanup_files(paths)

  ldscores <- Sys.getenv("siletti_superclusters_ldscores")

  pldsc <- cleansumstats_pldsc(paths,ldscores)



  captured <- c(
    slurm_header,
    "# 1) Cleansumstats pipeline, to QC and harmonize summary statistic",
    cleansumstats_job,
    "\n",
    "# Create directories, and calculate the number of significant SNPs",
    create_dirs,
    sig_snps,
    fix_pvals,
    "\n",
    "# 2) LDSC to munge and estimate heritability",
    ldsc_job,
    "\n",
    "# 3) Clumping to get number of loci",
    clump_job,
    "\n",
    "\n ",
    "# 5) Generate code to run SbayesR for polygenic scores ",
    # sbayes_job,
    # sbatch_sbayes,
    "\n",
    "# Clean up files that are no longer needed",
    cleanup
  )

  writeLines(captured, paths[["base_job"]])

  paths[["base_job"]]
}

#' Represent very small pvalues as minimum double number
#'
#' @param paths filepath for gwas file
#'
#' @return writes out the gwas file
#' @export
#'
#' @examples \dontrun{
#' correct_small_pval("/my_gwas/gwas.tsv.gz")
#' }
correct_small_pval <- function(paths) {

 first <- glue::glue("gunzip -c {paths[['clean']]} | awk -F '\\t' 'BEGIN{{OFS=\"\\t\"}}{{if (NR==1) {{for (i=1; i<=NF; i++) {{if ($i==\"P\") {{pcol=i; break}}}}}} else if ($pcol <  2.225074e-307) {{$pcol = 2.225074e-307}} {{print}}}}' | gzip > {paths[['temp']]}")
 second <- glue::glue("mv {paths[['temp']]} {paths[['clean']]}")
 c(first, second)

}

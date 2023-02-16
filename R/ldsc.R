run_ldsc <- function(paths) {

  add <- ldsc_col_checker(paths[["metafile"]])

  # Create the munge code
  base <- list(
    glue::glue("munge_sumstats.py "),
    glue::glue("--sumstats {paths[['clean']]} "),
    glue::glue("--out {paths[['ldsc_sumstats']]} "),
    glue::glue("--snp RSID "),
    glue::glue("--a1 EffectAllele "),
    glue::glue("--a2 OtherAllele "),
    glue::glue("--merge-alleles {Sys.getenv('ldsc_hm3')} "),
    glue::glue("--chunksize 500000 ")
  )


  # add in effect modifier and ncas+ncon or n arguments
  munge <-base %>%
    purrr::reduce(c) %>%
    c(., purrr::reduce(add,c)) %>%
    stringr::str_c(collapse = "")


  intercept <- glue::glue(
    "ldsc.py ",
    "--h2 {paths[['ldsc_sumstats']]}.sumstats.gz ",
    "--ref-ld-chr {Sys.getenv('ldsc_eur_wld')} ",
    "--w-ld-chr {Sys.getenv('ldsc_eur_wld')} ",
    "--out {paths[['ldsc_out']]}"

  )

  c(c("module unload python", "module load ldsc"),  munge, intercept)
}

ldsc_col_checker <- function(path) {
  metadata <- yaml::yaml.load_file(path)
  add <- list()

  # Check if to us B or Z for ldsc
  if("col_OR" %in% names(metadata) | "col_BETA" %in% names(metadata)) {
    add[[1]] <- glue::glue("--signed-sumstats B,0 ")
  } else {
    add[[1]] <- glue::glue("--signed-sumstats Z,0 ")
  }

  #  If possible use case-control sample size.
  if("col_CaseN" %in% names(metadata) & "col_ControlN" %in% names(metadata)){
    add <- c(add,
      glue::glue("--N-cas-col CaseN "),
      glue::glue("--N-con-col ControlN ")
      )
  } else {
    add <- c(add, glue::glue("--N-col N"))
  }

  add
}



cleansumstats_pldsc <- function(paths, ldscores) {

  all_celltypes <- fs::dir_ls(ldscores)


  job1 <- purrr::map2(
    all_celltypes,
    fs::path_file(all_celltypes),
    \(celltype, name) ldsc_partitioned(
      ld1 = celltype, outname = name, sumstats = paste0(paths[["ldsc_sumstats"]], ".sumstats.gz"),
      outdir = paths[["pldsc_siletti"]], base_ldscore = Sys.getenv("pldsc_base_ldscore"),
      weights = Sys.getenv("pldsc_weights"), freq = Sys.getenv("pldsc_freq"))
  )

  purrr::map_chr(job1, \(code) glue::glue("sbatch --time=00:30:00 --output={paths[['pldsc_siletti']]}/slurm-%j.out --mem=4gb --wrap='module unload python && module load ldsc && {code}'"))


}

ldsc_partitioned <- function(
    ld1,
    outname,
    sumstats,
    outdir,
    base_ldscore,
    weights,
    freq
)
  {


    glue::glue(
      "ldsc.py ",
      "--h2 {sumstats} ",
      "--ref-ld-chr {base_ldscore},{ld1}/baseline. ",
      "--w-ld-chr {weights} ",
      "--overlap-annot ",
      "--frqfile-chr {freq} ",
      "--print-coefficients ",
      "--out {outdir}/{outname}"
    )

  }






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

  c(c("module unload python", "module load ldsc"), "\n", munge, "\n", intercept)
}

ldsc_col_checker <- function(path) {
  metadata <- yaml::yaml.load_file(path)
  add <- list()

  # Check if to us B or Z for ldsc
  if(purrr::has_element(metadata, "BETA")) {
    add[[1]] <- glue::glue("--signed-sumstats B,0 ")
  } else {
    add[[1]] <- glue::glue("--signed-sumstats Z,0 ")
  }

  #  If possible use case-control sample size.
  if(purrr::has_element(metadata, "col_CaseN") & purrr::has_element(metadata, "col_CaseN")){
    add <- c(add,
      glue::glue("--N-cas-col CaseN "),
      glue::glue("--N-con-col ControlN ")
      )
  } else {
    add <- c(add, glue::glue("--N-col N"))
  }

  add
}





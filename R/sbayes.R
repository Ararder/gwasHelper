paths <- filepath_manager("scz2022_eur")
sbayes_cleansumstats <- function(paths) {
  # wrapper for generating commandline code
  code <- glue::glue("gwasHelper::prepare_ma('{paths[['clean']]}', '{paths[['sbayes_ma']]}')")
  cleaning <- glue::glue("Rscript -e {code}")

  sbayes_code <- gctb_sbayesr(
    infile = paths[["sbayes_ma"]],
    name = paths[['dataset_name']],
    out  = paths[['sbayes_dir']]
    )

  # make slurm header
  header <- slurm_header(24, 40, Sys.getenv("SLURM_JOB_ID"))

  # write slurm job

  writeLines(
    c(header, "\n", cleaning, "\n", sbayes_code),
    paths[["sbayes_slurm_path"]],
  )

  # system(glue("sbatch {paths[['sbayes_slurm_path']]"))
}

#' Title
#'
#' @param path filepath for sumstat
#' @param out filepath for .ma file
#' @param N Sample size - optional to use if missing N
#' @param freq_path path to SNP, A1, A2 freq if allele frequency is missing
#'
#' @return a character vector (commandline code)
#' @export
#'
#' @examples \dontrun{
#' prepare_ma(paths)
#' }
#'
prepare_ma <- function(path, out, N, freq_path=Sys.getenv("freq1000g")) {

  # read in file
  df <- data.table::fread(path)

  # apply changes if needed
  if(!missing(N)) {
    df <- sbayes_filters(df, N, freq_path)
  } else {
    df <- sbayes_filters(df, freq_path)
  }
  # Which allele frequencies to apply if missing?

  stopifnot("SbayesR requires SNP A1 A2 freq b se p N" = all(
    c("RSID", "EffectAllele","OtherAllele", "EAF", "B", "SE", "P", "N") %in% colnames(df)))

  # write out file

  dplyr::select(df,SNP = RSID, A1 = EffectAllele, A2 = OtherAllele, freq=EAF, b=B, se=SE, p=P, N=N) %>%
    readr::write_tsv(out)

}


sbayes_filters <- function(df, N, freq_path){
  # Apply filters on allele frequency and INFO if those columns exist

  if(!"EAF" %in% colnames(df)) {
    freq <- read_tsv(freq_path)
    df <- add_1000g_freq(df)
    print("Frequency missing. Adding 1000kg eur freq")
  }

  df <- filter(df, EAF >= 0.01 & EAF <= 0.99)

  if("INFO" %in% colnames(df)) {
    print("INFO detected. Filtering to INFO >= 0.9")
    df <- filter(df, INFO >= 0.9)
  }

  if(!missing(N)){
    df <- mutate(df, N = {{ N }})
  }

  if(all(c("ControlN", "CaseN") %in% colnames(df))) {
    df <- mutate(df, N = ControlN + CaseN)
  }

  df
}

add_1000g_freq <- function(df, freq) {

  matching <- inner_join(df,freq,by = c("RSID" = "SNP", "EffectAllele" = "A1", "OtherAllele"= "A2"))
  flipped <- inner_join(df,freq,by = c("RSID" = "SNP", "EffectAllele" = "A2","OtherAllele" = "A1")) %>%
    mutate(EAF = 1-EAF)

  bind_rows(matching,flipped)

}


gctb_sbayesr <- function(infile,
                             name,
                             pi="0.95,0.02,0.02,0.01",
                             gamma="0.0,0.01,0.1,1",
                             chain_length=10000,
                             burn_in=4000,
                             thin=10,
                             out_freq=10,
                             seed=2022,
                             out,
                             ldmatrix=Sys.getenv("GCTB_LDMATRIX"),
                             gctb=Sys.getenv("GCTB")
) {


  glue::glue(

    "{gctb} ",
    "--sbayes R ",
    "--mldm {ldmatrix} ",
    "--gwas-summary {infile} ",
    "--pi {pi} ",
    "--gamma {gamma} ",
    "--robust ",
    "--exclude-mhc ",
    "--impute-n ",
    "--seed {seed} ",
    "--chain-length {chain_length} ",
    "--burn-in {burn_in} ",
    "--thin {thin} ",
    "--out-freq {out_freq} ",
    "--no-mcmc-bin ",
    "--out {out}/{name}"
  )
}
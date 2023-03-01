utils::globalVariables(c("CHR", "magma_part1", "magma_geneset_sc", "FULL_NAME"))
#' Generate code to run magma gene assoc + gene_set assoc
#'
#' @param paths list of filepaths, defined using filepath_manager()
#' @param geneset_annotation_files character vector of filepaths to gene annotations
#'
#' @return a slurm job (as a character vector)
#' @export
#'
#' @examples \dontrun{
#' cleansumstats_magma(filepath_manager("test_gwas"))
#' }
#'
cleansumstats_magma <- function(paths,
                                geneset_annotation_files = c(
                                  Sys.getenv("siletti_magma_clusters"),
                                  Sys.getenv("siletti_magma_superclusters")
                                  )
                                ) {

  header <- slurm_header(24, 20, paths[["magma_dir"]])
  munge_magma <- glue::glue(
    "Rscript -e 'gwasHelper::magma_munge_and_setup(commandArgs(trailingOnly=TRUE)[2],commandArgs(trailingOnly=TRUE)[3])'",
    " --args {paths[['clean']]} {paths[['magma_dir']]}"
    )
  magma_gene <- magma_gene(dir = paths[["magma_dir"]])

  magma_genesets <- purrr::map(geneset_annotation_files, \(file) magma_geneset(
    dir = paths[["magma_dir"]],
    annotation_file = file,
    annotation_name = stringr::str_remove(fs::path_file(file), "_top10pct_genesets.txt")
    )) |>
    purrr::reduce(c)

  c(header, "module load magma", munge_magma, magma_gene, magma_genesets)
}




#' Create the nessecary field for MAGMA gene & gene-set analysis
#'
#' @param input_gwas filepath to input GWAS
#' @param output_dir output directory of results
#'
#' @return Returns nothing, but writes files to output dir when called
#' @export
#'
#' @examples \dontrun{
#' magma_munge_and_setup("my_gwas", "/gwas_project/magma/version2")
#' }
magma_munge_and_setup <- function(input_gwas, output_dir){
  df <- dplyr::tibble(data.table::fread(input_gwas))
  fs::dir_create(output_dir, recurse = TRUE)
  # Magma pipeline is configured for GRCh37.
  # Read in build 37 coordinates
  b37 <- dplyr::tibble(data.table::fread(fs::path(fs::path_dir(input_gwas), "cleaned_GRCh37.gz")))

  df <- df |>
    dplyr::select(-CHR, -POS) |>
    dplyr::inner_join(b37) |>
    magma_filters()
  print("Applied filters for MAGMA, and mapped to genome build 37")

  # check that all needed columns exist
  validate_columns("magma", df)
  print("All columns that are required are present")

  # write out snploc file
  dplyr::select(df, RSID, CHR, POS) %>%
    data.table::fwrite(fs::path(output_dir, "snplocs.tsv") ,sep ="\t", col.names = FALSE)

  # write out pval file
  dplyr::select(df,SNP=RSID, P, N) %>%
    data.table::fwrite(fs::path(output_dir, "pval.tsv") ,sep ="\t", col.names = FALSE)

  print(glue::glue("Wrote files snplocs.tsv and pval.tsv to {output_dir}"))

}




magma_filters <- function(df){
  # Apply filters on allele frequency and INFO if those columns exist
  if("EAF" %in% colnames(df)) {
    df <- dplyr::filter(df, EAF >= 0.01 & EAF <= 0.99)
  }

  if("INFO" %in% colnames(df)) {
    df <- dplyr::filter(df, INFO >= 0.9)
  }

  if(all(c("ControlN", "CaseN") %in% colnames(df))) {
    df <- dplyr::mutate(df, N = ControlN + CaseN)
  }

  df
}

magma_gene <- function(
  dir, window="35,10",
  gene_loc=Sys.getenv("gene_loc"),
  bfile = Sys.getenv("magma_bfile")
  ) {

  snp_to_genes <- list(
    glue::glue("magma "),
    glue::glue("--annotate window={window} "),
    glue::glue("--snp-loc {dir}/snplocs.tsv "),
    glue::glue("--gene-loc {gene_loc} "),
    glue::glue("--out {dir}/gwas")
  ) %>%
    purrr::reduce(c) %>%
    stringr::str_c(collapse="")

  est_gene_assoc <-
    list(
      glue::glue("magma "),
      glue::glue("--bfile {bfile} "),
      glue::glue("--pval {dir}/pval.tsv ncol=3 "),
      glue::glue("--gene-annot {dir}/gwas.genes.annot ") ,
      glue::glue("--out {dir}/gwas")
    ) %>%
    purrr::reduce(c) %>%
    stringr::str_c(collapse="")

  c(snp_to_genes, est_gene_assoc)

}





magma_geneset <- function(dir, annotation_file, annotation_name){
  glue::glue(
    "magma ",
    "--gene-results {dir}/gwas.genes.raw ",
    "--set-annot {annotation_file} ",
    "--out {dir}/{annotation_name}.gene_set "
  )

}

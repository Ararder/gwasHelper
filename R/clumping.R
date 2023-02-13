utils::globalVariables(c("POS", "tmp", "chr", "start", "end", "N", "P", "SNP"))

run_clumping <- function(infile, outdir) {
  stopifnot(!missing(outdir))
  # infile corresponds to a GWAS sumstats
  # outdir to

  plink <- clump_plink(infile, outdir)
  format_munge <- glue::glue("R -e 'gwasHelper::ranges_to_bed(commandArgs(trailingOnly = TRUE)[1],commandArgs(trailingOnly=TRUE)[2])'") %>%
    paste0(" --args", " ", fs::path(outdir, "clumps.clumped.ranges"), " ", fs::path(outdir, "clumps.bed"))

  bed_i <- fs::path(outdir, "clumps.bed")
  bed_out <- fs::path(outdir, "n_loci.bedtools")
  bed <- glue::glue("bedtools merge -d 50000 -i {bed_i} -c 4,5 -o sum,min > {bed_out}")

  c("module load bedtools", "module load plink", plink, "\n", format_munge, "\n", bed)

}



clump_plink <- function(
  sumstat,
  outdir,
  p1 = "5e-08",
  p2 = "5e-06",
  r2 = 0.1,
  kb = 3000,
  snp_field = "RSID" ,
  p_field = "P",
  ref = Sys.getenv("genome1000"),
  range = Sys.getenv("ref_gene_list")
) {
  glue::glue(
    "plink --bfile {ref} ",
    "--clump {sumstat} ",
    "--out {outdir}/clumps ",
    "--clump-p1 {p1} ",
    "--clump-p2 {p2} ",
    "--clump-r2 {r2} ",
    "--clump-kb {kb} ",
    "--clump-snp-field {snp_field} ",
    "--clump-field {p_field} ",
    "--clump-range {range} "
  )

}


#' Conver the clumping.ranges format from plink --clump to bed compatible format
#'
#' @param infile a clumping.ranges file from plink --clump
#' @param out filename of output
#'
#' @return a character vector of commandline code
#' @export
#'
#' @examples \dontrun{
#' ranges_to_bed("analysis/plink/clumping/scz2022.clumping.ranges", "scz2022_clumping/ranges.bed")
#' }
ranges_to_bed <- function(infile, out){

data.table::fread(infile) %>%
  dplyr::mutate(
    chr = stringr::word(POS, 1, sep = stringr::fixed(":")),
    tmp = stringr::word(POS, 2, sep = stringr::fixed(":")),
    start = as.integer(stringr::word(tmp, 1, sep = stringr::fixed(".."))),
    end = as.integer(stringr::word(tmp, 2, sep = stringr::fixed("..")))
  ) %>%
  dplyr::arrange(chr, start, end) %>%
  dplyr::select(chr, start, end, N, P, SNP) %>%
  readr::write_tsv(out, col_names = FALSE)
}




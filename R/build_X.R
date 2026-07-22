#' Build a genotype matrix from chromosome-specific PLINK files
#'
#' Extracts requested SNPs from chromosome-specific PLINK BED files, merges
#' the extracted files, fills missing genotypes with the A2 allele, and returns
#' a dense genotype matrix.
#'
#' @param SNPInfo Data frame containing columns `CHR` and `SNP`.
#' @param bfile_pattern PLINK file prefix containing the placeholder
#'   \code{\{CHR\}}.
#' @param keep_file Optional PLINK keep file.
#' @param plink Path to the PLINK executable.
#' @param temp_root Directory used for temporary PLINK files.
#' @param clumping Whether to remove highly correlated SNPs with PLINK
#'   clumping. All SNPs receive the same p-value, and the clump radius is fixed
#'   at 1000 kb, so this performs LD filtering only.
#' @param r2 PLINK clumping r-squared threshold. Used only when `clumping` is
#'   `TRUE`.
#' @param verbose Whether to display PLINK output.
#'
#' @return A list containing `X`, the SNP information retained by PLINK,
#'   `FAM`, `SNP_missing`, and `SNP_clumped`. `SNP_missing` contains requested
#'   SNPs absent from the source BED files, while `SNP_clumped` contains SNPs
#'   removed by LD clumping. Rows of `X` correspond to rows of `FAM`.
#' @export
build_X <- function(SNPInfo, bfile_pattern, keep_file = NULL,
                    plink = "./plink", temp_root = tempdir(),
                    clumping = FALSE, r2 = 0.81, verbose = FALSE) {
  if (!is.data.frame(SNPInfo) ||
      !all(c("CHR", "SNP") %in% names(SNPInfo))) {
    stop("SNPInfo must contain columns named CHR and SNP.")
  }
  if (!grepl("{CHR}", bfile_pattern, fixed = TRUE)) {
    stop("bfile_pattern must contain {CHR}.")
  }
  if (!is.logical(clumping) || length(clumping) != 1L || is.na(clumping)) {
    stop("clumping must be TRUE or FALSE.")
  }
  if (clumping && (!is.numeric(r2) || length(r2) != 1L ||
                   !is.finite(r2) || r2 <= 0 || r2 > 1)) {
    stop("r2 must be a single number greater than 0 and no greater than 1.")
  }

  selected <- unique(data.frame(
    CHR = as.character(SNPInfo$CHR),
    SNP = as.character(SNPInfo$SNP),
    stringsAsFactors = FALSE
  ))
  if (!nrow(selected) || anyNA(selected) ||
      any(!nzchar(selected$CHR)) || any(!nzchar(selected$SNP))) {
    stop("SNPInfo contains no usable chromosome and SNP records.")
  }

  bfile_pattern <- path.expand(bfile_pattern)
  keep_file <- if (is.null(keep_file)) NULL else path.expand(keep_file)
  plink <- path.expand(plink)
  temp_root <- path.expand(temp_root)

  tmp <- tempfile(pattern = "build_X_", tmpdir = temp_root)
  dir.create(tmp, recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)

  run_plink <- function(args) {
    output <- if (verbose) "" else FALSE
    status <- system2(plink, args = args, stdout = output, stderr = output)
    if (status != 0L) {
      stop("PLINK exited with status ", status, ".", call. = FALSE)
    }
  }

  chr_levels <- unique(selected$CHR)
  chr_prefix <- character(length(chr_levels))
  for (i in seq_along(chr_levels)) {
    chr <- chr_levels[i]
    chr_tag <- gsub("[^A-Za-z0-9._-]", "_", chr)
    extract_file <- file.path(tmp, paste0("snps_", chr_tag, ".txt"))
    data.table::fwrite(
      selected[selected$CHR == chr, "SNP", drop = FALSE],
      extract_file,
      col.names = FALSE
    )

    chr_prefix[i] <- file.path(tmp, paste0("chr_", chr_tag))
    args <- c(
      "--bfile", shQuote(gsub("{CHR}", chr, bfile_pattern, fixed = TRUE)),
      "--extract", shQuote(extract_file)
    )
    if (!is.null(keep_file)) {
      args <- c(args, "--keep", shQuote(keep_file))
    }
    args <- c(args, "--make-bed", "--out", shQuote(chr_prefix[i]))
    run_plink(args)
  }

  merged_missing <- chr_prefix[1]
  if (length(chr_prefix) > 1L) {
    merge_file <- file.path(tmp, "merge_list.txt")
    writeLines(chr_prefix[-1], merge_file)
    merged_missing <- file.path(tmp, "merged_missing")
    run_plink(c(
      "--bfile", shQuote(chr_prefix[1]),
      "--merge-list", shQuote(merge_file),
      "--make-bed", "--out", shQuote(merged_missing)
    ))
  }

  bim_available <- data.table::fread(
    paste0(merged_missing, ".bim"),
    header = FALSE
  )
  SNP_available <- as.character(bim_available[[2]])
  SNP_missing <- setdiff(selected$SNP, SNP_available)
  SNP_clumped <- character()

  if (clumping) {
    clump_input <- file.path(tmp, "clump_input.txt")
    data.table::fwrite(
      data.frame(SNP = SNP_available, P = rep(1e-6, length(SNP_available))),
      clump_input,
      sep = "\t"
    )

    clump_prefix <- file.path(tmp, "clump")
    run_plink(c(
      "--bfile", shQuote(merged_missing),
      "--clump", shQuote(clump_input),
      "--clump-snp-field", "SNP",
      "--clump-field", "P",
      "--clump-kb", "1000",
      "--clump-p1", "1",
      "--clump-p2", "1",
      "--clump-r2", format(r2, scientific = FALSE, trim = TRUE),
      "--out", shQuote(clump_prefix)
    ))

    clumped <- data.table::fread(paste0(clump_prefix, ".clumped"))
    SNP_keep <- unique(as.character(clumped$SNP))
    SNP_clumped <- setdiff(SNP_available, SNP_keep)
    clump_extract <- file.path(tmp, "clump_extract.txt")
    data.table::fwrite(
      data.frame(SNP = SNP_keep),
      clump_extract,
      col.names = FALSE
    )

    clumped_bed <- file.path(tmp, "clumped")
    run_plink(c(
      "--bfile", shQuote(merged_missing),
      "--extract", shQuote(clump_extract),
      "--make-bed", "--out", shQuote(clumped_bed)
    ))
    merged_missing <- clumped_bed
  }

  merged <- file.path(tmp, "merged")
  run_plink(c(
    "--bfile", shQuote(merged_missing),
    "--fill-missing-a2", "--make-bed", "--out", shQuote(merged)
  ))

  bim <- data.table::fread(paste0(merged, ".bim"), header = FALSE)
  fam <- data.table::fread(
    paste0(merged, ".fam"),
    header = FALSE,
    colClasses = "character"
  )
  SNPInfo_out <- data.frame(
    CHR = as.character(bim[[1]]),
    SNP = as.character(bim[[2]]),
    BP = as.integer(bim[[4]]),
    A1 = as.character(bim[[5]]),
    A2 = as.character(bim[[6]]),
    stringsAsFactors = FALSE
  )
  FAM <- data.frame(
    FID = fam[[1]],
    IID = fam[[2]],
    PAT = fam[[3]],
    MAT = fam[[4]],
    SEX = fam[[5]],
    PHENO = fam[[6]],
    stringsAsFactors = FALSE
  )

  X <- as.matrix(BEDMatrix::BEDMatrix(paste0(merged, ".bed")))
  colnames(X) <- SNPInfo_out$SNP
  rownames(X) <- NULL

  list(
    X = X,
    SNPInfo = SNPInfo_out,
    FAM = FAM,
    SNP_missing = SNP_missing,
    SNP_clumped = SNP_clumped
  )
}

#' EURblock LD Reference Panel
#'
#' LD reference panel block derived from UK Biobank European population data
#' for SBayesRC method.
#'
#' @format A list with 6 elements:
#' \describe{
#'   \item{m}{Integer. Number of SNPs in this LD block (7346).}
#'   \item{k}{Integer. Number of components retained in low-rank approximation (2073).}
#'   \item{sumLambda}{Numeric vector. Cumulative sum of eigenvalues (length 7346).}
#'   \item{thresh}{Numeric. Threshold for cumulative proportion of eigenvalues (0.995).}
#'   \item{lambda}{Numeric vector. Eigenvalues from eigendecomposition (length 2073).}
#'   \item{U}{Numeric matrix. Eigenvectors matrix (7346 x 2073).}
#' }
#'
#' @source \url{https://github.com/zhilizheng/SBayesRC}
#'
#' @examples
#' \dontrun{
#' # Load the data
#' data(EURblock)
#'
#' # Access components
#' n_snps <- EURblock$m
#' n_components <- EURblock$k
#' }
"EURblock"

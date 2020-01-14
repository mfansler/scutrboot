#' two_sample_test_lui
#'
#' @param object
#' @param idx.A
#' @param idx.B
#' @param idx.SU
#' @param idx.LU
#' @param replicates
#' @param pseudocount
#'
#' @return
#' @export
#'
#' @examples
two_sample_test_lui <- function (object, idx.A, idx.B, idx.SU, idx.LU, replicates=2000, pseudocount=1) {
  cts <- to_counts(object)

  idx.AB <- c(idx.A, idx.B)

  n.A <- length(idx.A)
  n.B <- length(idx.B)
  n.AB <- n.A + n.B

  LUI.A.hat <- mean_lu_usage(cts[, idx.A], idx.SU, idx.LU)
  LUI.B.hat <- mean_lu_usage(cts[, idx.B], idx.SU, idx.LU)
  LUI.diff.hat <- LUI.B.hat - LUI.A.hat

  # bootstrap matrix for pseudo group A
  M.bs.A <- Matrix::Matrix(rmultinom(replicates, n.A, rep(1, n.AB)), sparse=TRUE)

  SU.bs.A  <- cts[idx.SU, idx.AB] %*% M.bs.A
  LU.bs.A  <- cts[idx.LU, idx.AB] %*% M.bs.A
  LUI.bs.A <- LU.bs.A / (SU.bs.A + LU.bs.A)

  # bootstrap matrix for pseudo group B
  M.bs.B <- Matrix::Matrix(rmultinom(replicates, n.B, rep(1, n.AB)), sparse=TRUE)

  SU.bs.B  <- cts[idx.SU, idx.AB] %*% M.bs.B
  LU.bs.B  <- cts[idx.LU, idx.AB] %*% M.bs.B
  LUI.bs.B <- LU.bs.B / (SU.bs.B + LU.bs.B)

  LUI.bs.pvals <- Matrix::rowMeans((abs(LUI.bs.B - LUI.bs.A) - abs(LUI.diff.hat)) >= 0)

  # adjust by pseudocount
  LUI.bs.pvals <- (LUI.bs.pvals*replicates + pseudocount) / (replicates + pseudocount)

  data.frame(A.hat=LUI.A.hat, B.hat=LUI.B.hat, pvals=LUI.bs.pvals)
}

#' two_sample_test_gene
#'
#' @param object
#' @param idx.A
#' @param idx.B
#' @param replicates
#' @param pseudocount
#'
#' @return
#' @export
#'
#' @examples
two_sample_test_gene <- function (object, idx.A, idx.B, replicates=2000, pseudocount=1) {
  cts <- to_counts(object)

  idx.AB <- c(idx.A, idx.B)

  n.A <- length(idx.A)
  n.B <- length(idx.B)
  n.AB <- n.A + n.B

  gene.A.hat <- mean_counts_per_cell(cts[, idx.A])
  gene.B.hat <- mean_counts_per_cell(cts[, idx.B])
  gene.diff.hat <- gene.B.hat - gene.A.hat

  # bootstrap matrix for pseudo group A
  M.bs.A <- Matrix::Matrix(rmultinom(replicates, n.A, rep(1, n.AB)), sparse=TRUE)

  gene.bs.A  <- cts[, idx.AB] %*% M.bs.A / n.A

  # bootstrap matrix for pseudo group B
  M.bs.B <- Matrix::Matrix(rmultinom(replicates, n.B, rep(1, n.AB)), sparse=TRUE)

  gene.bs.B  <- cts[, idx.AB] %*% M.bs.B / n.B

  gene.bs.pvals <- Matrix::rowMeans((abs(gene.bs.B - gene.bs.A) - abs(gene.diff.hat)) >= 0)

  # adjust by pseudocount
  gene.bs.pvals <- (gene.bs.pvals*replicates + pseudocount) / (replicates + pseudocount)

  data.frame(A.hat=gene.A.hat, B.hat=gene.B.hat, pvals=gene.bs.pvals)
}

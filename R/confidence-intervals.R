#' lui_ci
#'
#' @param object
#' @param idx.SU
#' @param idx.LU
#' @param replicates
#' @param probs
#' @param na.rm
#'
#' @return
#'
#' @import Matrix
#' @importFrom sparseMatrixStats rowQuantiles
#' @export
lui_ci <- function(object, idx.SU, idx.LU, replicates = 2000, probs = c(0.005, 0.025, 0.25, 0.5, 0.75, 0.975, 0.995), na.rm = TRUE) {
    cts <- to_counts(object)

    n.cells <- ncol(cts)

    # bootstrap matrix
    M.bs <- Matrix(rmultinom(replicates, n.cells, rep(1, n.cells)), sparse = TRUE)

    SU.bs <- cts[idx.SU, ] %*% M.bs
    LU.bs <- cts[idx.LU, ] %*% M.bs
    LUI.bs <- LU.bs / (SU.bs + LU.bs)

    qs <- rowQuantiles(LUI.bs, probs = probs, na.rm = na.rm)

    qs
}

#' gene_ci
#'
#' @param object
#' @param replicates
#' @param probs
#' @param na.rm
#'
#' @return
#'
#' @import Matrix
#' @importFrom sparseMatrixStats rowQuantiles
#' @export
gene_ci <- function(object, replicates = 2000, probs = c(0.005, 0.025, 0.25, 0.5, 0.75, 0.975, 0.995), na.rm = TRUE) {
    cts <- to_counts(object)

    n.cells <- ncol(cts)

    # bootstrap matrix
    M.bs <- Matrix(rmultinom(replicates, n.cells, rep(1, n.cells)), sparse = TRUE)

    gene.bs <- as.matrix(cts %*% M.bs / n.cells)

    qs <- rowQuantiles(gene.bs, probs = probs, na.rm = na.rm)
    qs
}

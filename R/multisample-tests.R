## computes the sum of squared residuals of a matrix row
.compute_rss <- function(lui, na.rm = FALSE) {
    if (length(dim(lui)) < 2) {
        sum((lui - mean(lui, na.rm = na.rm))^2, na.rm = na.rm)
    } else {
        rowSums((lui - rowMeans(lui, na.rm = na.rm))^2, na.rm = na.rm)
    }
}

## internal serialization method
## converts a matrix row into a string representation of column combinations
.row_to_idxstr <- function(r) {
    paste(which(r), collapse = ",")
}

## internal deserialization method
## converts string representation to column index vector
.idxstr_to_idx <- function(idxstr) {
    as.integer(unlist(strsplit(idxstr, split = ",", fixed = TRUE)))
}

.chunk <- function(x, maxChunkSize) {
    n_chunks <- ceiling(length(x) / maxChunkSize)
    if (n_chunks > 1) {
          split(x, cut(seq_along(x), n_chunks, labels = FALSE))
      } else {
          x
      }
}

## Credit: https://stackoverflow.com/a/49252734/570918
.flatten <- function(x, use.names = TRUE, classes = "ANY") {
    #' Source taken from rlist::list.flatten
    len <- sum(rapply(x, function(x) 1L, classes = classes))
    y <- vector("list", len)
    i <- 0L
    items <- rapply(x, function(x) {
        i <<- i + 1L
        y[[i]] <<- x
        TRUE
    }, classes = classes)
    if (use.names && !is.null(nm <- names(items))) {
          names(y) <- nm
      }
    y
}

.chunk_oversized <- function(groups_to_genes_dict, maxChunkSize) {
    ## break up chunks that are larger than maxChunkSize
    chunked_dict <- lapply(groups_to_genes_dict, .chunk, maxChunkSize = maxChunkSize)

    ## flatten, keeping names
    chunked_dict <- .flatten(chunked_dict)

    ## restore names to original name set
    names(chunked_dict) <- gsub(".[0-9]+$", "", names(chunked_dict))

    chunked_dict
}

#' Perform multisample bootstrap statistical test on variance in LUI estimates.
#'
#' @param sce A \code{SingleCellExperiment} object with 3' UTR isoform counts.
#' @param groupKey A string for the \strong{group} key corresponding to a
#'   column in \code{colData(sce)}.
#' @param isoformKey A string for the \strong{isoform} key corresponding to a
#'   column in \code{rowData(sce)}.
#' @param geneKey A string for the \strong{gene} key corresponding to a column
#'   in \code{rowData(sce)}.
#' @param minCells An integer threshold for minimum number of cells in a group
#'   that express a gene. Default is 50 cells.
#' @param minCoexpressed An integer threshold for minimum number of groups
#'   that satisfied the \code{minCells} threshold. Must be at least 2 (default).
#' @param nPermutations An integer number of permutations to perform.
#' @param enforceMinCells A boolean specifying whether to discard permutations
#'   not satisfying the \code{minCells} threshold in all groups being tested.
#' @param maxChunkSize An integer specifying the maximum number of genes to
#'   process in parallel. Smaller chunks can help to better distribute the test
#'   across multiple cores.
#' @param assayName the name of the assay with UMI counts
#'
#' @return A data.frame summarizing the test results for all genes that could be
#'   tested given the specified \code{minCells} and \code{minCoexpressed}
#'   values. The rownames and first column correspond to \strong{gene} names.
#'   The subsequent columns are named after the \strong{group} levels, and
#'   contain the point estimates for the LU index values. These are followed by:
#'    - \code{n_groups_minCells}, the count of groups that satisfied the
#'   \code{minCells} threshold;
#'    - \code{max_lui_diff}, the maximum difference in points estimates;
#'    - \code{p_val}
#'    - \code{q_val}, a Benjamini-Hochberg-adjusted p-value
#'
#'
#' @import SingleCellExperiment
#' @importFrom Matrix fac2sparse
#' @importFrom sparseMatrixStats rowMaxs rowMins
#' @importFrom BiocParallel bpmapply
#' @export
multisample_test_lui <- function(sce, groupKey = "cell_type", isoformKey = "isoform", geneKey = "gene",
    minCells = 50, minCoexpressed = 2, nPermutations = 10000, enforceMinCells = TRUE,
    maxChunkSize = 30, assayName="counts") {
    ## DESIGN MATRICES
    ## ===============

    ## [cells x groups]
    M.groups <- t(fac2sparse(colData(sce)[[groupKey]]))
    rownames(M.groups) <- colnames(sce)

    ## [genes x transcripts]
    M.genes <- fac2sparse(rowData(sce)[[geneKey]])
    colnames(M.genes) <- rownames(sce)

    ## [genes x transcripts]
    M.LU <- M.genes %*% Diagonal(ncol(M.genes), rowData(sce)[[isoformKey]] == "LU")

    ## COMPUTE POINT ESTIMATES
    ## =======================

    ## extract raw counts
    cts.txs <- assay(sce, assayName)

    ## aggregate to group counts [transcripts x groups]
    cts.groups <- cts.txs %*% M.groups

    ## compute LUI point estimates for groups [genes x groups]
    lui.hat <- (M.LU %*% cts.groups) / (M.genes %*% cts.groups)

    ## compute number of cells expressing [genes x groups]
    expr.gene.groups <- ((M.genes %*% cts.txs) > 0) %*% M.groups

    ## mask LUI values for groups below threshold
    lui.hat[which(expr.gene.groups < minCells)] <- NA

    ## determine testable genes [coexpr_genes x 1]
    idx.coexpr <- which(rowSums(!is.na(lui.hat)) >= minCoexpressed)

    ## compute sum of squared residuals for testable genes [coexpr_genes x 1]
    rss.hat <- .compute_rss(lui.hat[idx.coexpr, ], na.rm = TRUE)

    ## extract unique combinations of groups as string serializations
    genes.groups <- apply(!is.na(lui.hat), 1, .row_to_idxstr)

    ## construct map {group_combination -> coexpr_genes}
    coexpr.dict.genes <- split(names(genes.groups), genes.groups)

    ## chunk dictionary for improved parallelization
    coexpr.dict.genes <- .chunk_oversized(coexpr.dict.genes, maxChunkSize)

    ## RSS PERMUTATIONS
    ## ================
    compute_rss_permutations <- function(genes, groups_idxstr) {
        ## deserialize group_combinations index
        groups <- .idxstr_to_idx(groups_idxstr)

        ## base return to catch unexpressed genes
        if (length(groups) < minCoexpressed) { # NOTE: this is redundant!
            return(setNames(rep(NA_real_, length(genes)), genes))
        }

        ## extract just the cells for group combination [cells x 1]
        cells <- which(rowSums(M.groups[, groups]) > 0)

        ## slice out design matrix
        M.subgroups <- M.groups[cells, groups]

        ## generate all permutations [cells x permutations]
        M.perms <- replicate(nPermutations, sample(seq_along(cells)))

        ## get indices of transcripts for the genes
        txs <- which((if (length(genes) > 1) colSums(M.genes[genes, ]) else M.genes[genes, ]) > 0)

        ## extract only the transcript counts for the cells being permuted
        cts.txs <- assay(sce[txs, cells], assayName)

        ## for each permutation (column in M.perms)
        matrix(apply(M.perms, MARGIN = 2, function(idx.perm) {
            ## aggregate permuted transcript counts by group
            cts.groups <- cts.txs %*% M.subgroups[idx.perm, ]

            ## compute group-gene LUIs
            lui.groups <- (M.LU[genes, txs] %*% cts.groups) / (M.genes[genes, txs] %*% cts.groups)
            rownames(lui.groups) <- genes

            ## if enforcing minimum number of cells per group in permutations
            if (enforceMinCells) {
                ## compute number of non-zero cells
                nnz_cells <- ((M.genes[genes, txs] %*% cts.txs) > 0) %*% M.subgroups[idx.perm, ]

                # mask any entry where a group had less than minimum number of cells expressing
                lui.groups[which(nnz_cells < minCells)] <- NA
            }

            ## return the sum of squared residuals
            .compute_rss(lui.groups)
        }), ncol = nPermutations, byrow = FALSE, dimnames = list(genes))
    }

    ## apply the above function to each group combination (returns list of matrices)
    ## then bind the matrices
    rss.perms <- do.call(
        rbind,
        bpmapply(compute_rss_permutations,
            genes = coexpr.dict.genes,
            groups_idxstr = names(coexpr.dict.genes),
            SIMPLIFY = FALSE, USE.NAMES = FALSE
        )
    )

    ## validate testing
    idx.failure <- names(rss.hat)[!(names(rss.hat) %in% rownames(rss.perms))]
    if (length(idx.failure) > 0) {
        warning(
            "\nThe following genes failed to generate valid bootstraps:\n",
            paste(idx.failure, collapse = "\n")
        )
    }

    idx.success <- names(rss.hat)[names(rss.hat) %in% rownames(rss.perms)]

    ## compute p-values
    ## NB: Includes pseudocount (consider the actual observation as a random one).
    ## This corresponds to a conservative adjustment to the p-value which
    ## transforms this to an estimator on the upper bound of the true p-value.
    p_vals <- (rowSums(rss.perms[idx.success, ] >= rss.hat[idx.success], na.rm = TRUE) + 1) / (rowSums(!is.na(rss.perms[idx.success, ])) + 1)

    ## compute q-values
    q_vals <- p.adjust(p_vals, "fdr")

    ## more useful at this point as a standard base::matrix object
    lui.hat <- as.matrix(lui.hat[idx.coexpr, ][idx.success, ])

    ## bind together a data.frame
    cbind(
        gene = rownames(lui.hat), as.data.frame(lui.hat),
        n_groups_minCells = rowSums(!is.na(lui.hat)),
        max_lui_diff = rowMaxs(lui.hat, na.rm = TRUE) - rowMins(lui.hat, na.rm = TRUE),
        p_val = p_vals, q_val = q_vals
    )
}

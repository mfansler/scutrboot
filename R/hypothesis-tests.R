#' Two Sample Test
#'
#' @param sce a \code{SingleCellExperiment} object
#' @param sampleKey a string specifying the factor column in \code{colData(sce)}
#'     that codes for samples
#' @param sample0 a string specifying the level for the baseline sample
#' @param sample1 a string specifying the level for the contrasting sample
#' @param statistic either "WD" (Wasserstein Distance) or "UI" (Usage Index). If
#'     Usage Index is specified, then the \code{featureIndex} argument must also
#'     be provided.
#' @param nBootstraps number of bootstraps samples to use estimate pvalues
#' @param geneKey a string specifying the factor column in \code{rowData(sce)}
#'     that specifies the gene with which each transcript is associated
#' @param featureIndex a string specifying a logical column in \code{rowData(sce)}
#'     that specifies which transcripts are indexed. The logical should be TRUE
#'     for at least one isoform per gene, otherwise the index will constant.
#' @param featureExclude a string specifying a logical column in \code{rowData(sce)}
#'     of features to exclude. All isoforms with a TRUE value in the column will
#'     be removed.
#' @param minCellsPerGene a numeric specifying the minimum number of cells
#'     that express a gene to consider the gene testable.
#' @param assay the assay in the \code{SingleCellExperiment} to use
#'
#' @return a \code{DataFrame} with columns
#'     - \code{gene}
#'     - \code{stat}
#'     - \code{pval}
#'
#'
#' @importFrom Matrix fac2sparse Diagonal drop0
#' @importFrom sparseMatrixStats rowAlls colAnys
#' @import SingleCellExperiment
#' @export
testTwoSample <- function (sce, sampleKey, sample0, sample1,
                           statistic=c("WD", "UI"), nBootstraps=10000,
                           geneKey="gene_id", featureIndex=NULL,
                           featureExclude=NULL, minCellsPerGene=50,
                           assayName="counts", pseudocount=1) {
    if (statistic == 'UI' && is.null(featureIndex)) {
        stop("The 'UI' (Usage Index) statistic requires a `featureIndex` to be provided!")
    } else if (statistic == 'WD' && !is.null(featureIndex)) {
        warning("The 'WD' (Wasserstein Distance) statistic ignores the `featureIndex` argument.")
    }

    ## filter cells
    sce <- sce[,colData(sce)[[sampleKey]] %in% c(sample0, sample1)]

    ## filter excluded transcripts
    if (!is.null(featureExclude)) {
        if (is.logical(featureExclude)) {
            sce <- sce[!featureExclude,]
        } else if (featureExclude %in% names(rowData(sce))) {
            sce <- sce[!rowData(sce)[[featureExclude]],]
        } else {
            warning("The `featureExclude` argument must either be a logical or",
                    " a column name of `rowData(sce)`! Ignoring argument.")
        }
    }

    ## filter unexpressed transcripts
    sce <- sce[rowSums(assay(sce, assayName)) > 0,]

    ## filter untestable transcripts
    cts_tx_cell <- assay(sce, assayName)

    M_gene_tx <- fac2sparse(rowData(sce)[[geneKey]])
    M_cell_sample <- t(fac2sparse(colData(sce)[[sampleKey]]))[,c(sample0, sample1)]

    ncells_gene_sample <- ((M_gene_tx %*% cts_tx_cell) > 0) %*% M_cell_sample

    idxExpressedGenes <- rowAlls(drop0(ncells_gene_sample >= minCellsPerGene))
    idxMultiTxGenes <- rowSums(M_gene_tx) > 1
    idxTestableTxs <- colAnys(M_gene_tx[idxExpressedGenes & idxMultiTxGenes,])

    sce <- sce[idxTestableTxs,]

    rm(list=c("ncells_gene_sample", "idxExpressedGenes", "idxMultiTxGenes",
              "idxTestableTxs"))

    ## compute observed statistics
    cts_tx_cell <- assay(sce, assayName)

    M_gene_tx <- fac2sparse(rowData(sce)[[geneKey]]) ## updated

    cts_tx_sample <- cts_tx_cell %*% M_cell_sample
    cts_gene_sample <- M_gene_tx %*% cts_tx_sample
    usage_tx_sample <- (t(M_gene_tx) %*% (1/cts_gene_sample)) * cts_tx_sample

    dUsage_tx <- usage_tx_sample[,sample1] - usage_tx_sample[,sample0]
    if (statistic == "WD") {
        stat_gene <- M_gene_tx %*% abs(dUsage_tx) / 2.0
    } else {
        D_index <- Diagonal(nrow(sce), rowData(sce)[[featureIndex]])
        stat_gene <- M_gene_tx %*% D_index %*% dUsage_tx
    }

    ## bootstrap statistics
    ncells_sample <- colSums(M_cell_sample)
    M_cell_bs0 <- Matrix(sparse=TRUE,
                         data=rmultinom(nBootstraps,
                                        ncells_sample[[sample0]],
                                        rep(1, sum(ncells_sample))))
    cts_tx_bs0 <- cts_tx_cell %*% M_cell_bs0
    cts_gene_bs0 <- M_gene_tx %*% cts_tx_bs0
    usage_tx_bs0 <- (t(M_gene_tx) %*% (1/cts_gene_bs0)) * cts_tx_bs0
    rm(M_cell_bs0, cts_tx_bs0, cts_gene_bs0); gc()

    M_cell_bs1 <- Matrix(sparse=TRUE,
                         data=rmultinom(nBootstraps,
                                        ncells_sample[[sample1]],
                                        rep(1, sum(ncells_sample))))
    cts_tx_bs1 <- cts_tx_cell %*% M_cell_bs1
    cts_gene_bs1 <- M_gene_tx %*% cts_tx_bs1
    usage_tx_bs1 <- (t(M_gene_tx) %*% (1/cts_gene_bs1)) * cts_tx_bs1
    rm(M_cell_bs1, cts_tx_bs1, cts_gene_bs1, cts_tx_cell); gc()

    ## Note: switching order of Kronecker product effectively computes all
    ##       combinations of columns
    #K_repeat <- matrix(rep(1, nBootstraps), nrow=1)
    #dUsage_tx_bs <- (K_repeat %x% usage_tx_bs1) - (usage_tx_bs0 %x% K_repeat)
    dUsage_tx_bs <- usage_tx_bs1 - usage_tx_bs0
    rm(usage_tx_bs1, usage_tx_bs0); gc()

    if (statistic == "WD") {
        stat_gene_bs <- M_gene_tx %*% abs(dUsage_tx_bs) / 2.0
    } else {
        stat_gene_bs <- M_gene_tx %*% D_index %*% dUsage_tx_bs
    }
    rm(dUsage_tx_bs); gc()

    pval_gene <- rowMeans(abs(stat_gene_bs) - abs(stat_gene) >= 0, na.rm=TRUE)
    nbs_true <- rowSums(!is.na(stat_gene_bs))

    ## adjust by pseudocount
    ##pval_gene <- (pval_gene * nBootstraps^2 + pseudocount) / (nBootstraps^2 + pseudocount)
    pval_gene <- (pval_gene * nbs_true + pseudocount) / (nbs_true + pseudocount)

    ## format dataframe
    DataFrame(gene=rownames(M_gene_tx), stat=stat_gene[,1],
              pval=pval_gene, bootstraps=nbs_true)
}

mean_counts_per_cell <- function(object, nonzero = FALSE, na.rm = TRUE) {
    cts <- to_counts(object)

    if (nonzero) {
        Matrix::rowSums(cts, na.rm = na.rm) / Matrix::rowSums(cts > 0, na.rm = na.rm)
    } else {
        Matrix::rowMeans(cts, na.rm = na.rm)
    }
}

mean_lu_usage <- function(object, idx.SU, idx.LU, na.rm = TRUE) {
    cts <- to_counts(object)
    SU <- Matrix::rowSums(cts[idx.SU, ], na.rm = na.rm)
    LU <- Matrix::rowSums(cts[idx.LU, ], na.rm = na.rm)
    LU / (SU + LU)
}

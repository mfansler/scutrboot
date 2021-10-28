to_counts <- function(object) {
    # verify class
    if (is(object, "SingleCellExperiment")) {
        counts(object)
    } else if (is(object, "Matrix") | is(object, "matrix")) {
        object
    } else {
        stop(sprintf("Object must be of class 'Matrix' or 'SingleCellExperiment': found %s", class(object)))
    }
}

sum_rows_by <- function(sce, by) {
    if (!(by %in% colnames(rowData(sce)))) {
        stop(sprintf("Column '%s' was not found in `rowData`. Available columns are:\n%s", by, colnames(rowData(sce))))
    }
    M.factors <- fac2sparse(rowData(sce)[, by])
    M.factors %*% counts(sce)
}

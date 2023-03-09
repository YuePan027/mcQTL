#' Estimation of cellular expression in high-throughput data from heterogeneous tissues
#'
#' This function returns a `SummarizedExperiment` object including deconvoluted cellular sample-wise expression.
#'
#' This is a function developed to implement deconvolution for cell-type-specific expression per bulk sample using `TCA`.
#'
#' @param se A `SummarizedExperiment` object with bulk protein expression data frame contained in `counts` slot, and
#' cell-type proportion estimates for each sample contained as an element (`prop`) in `metadata` slot.
#'
#' @return A `SummarizedExperiment`. The results after `TCA` deconvolution will be stored as an element (`TCA_deconv`) in `metadata` slot.
#' It is a list with the length of the number of cell types (as in `prop` in `metadata` slot).
#' Each element stores the deconvoluted protein expression matrix per bulk sample within that cell type.
#'
#' @import SummarizedExperiment
#'
#' @export
#'
#'
TCA_deconv <- function(se){
  res_tca <- TCA::tca(X = as.matrix(assay(se)),
                      W = se@metadata$prop,
                      refit_W = FALSE,
                      verbose = FALSE)
  res_full <- TCA::tensor(tca.mdl = res_tca, X = as.matrix(assay(se)))
  names(res_full) <- colnames(se@metadata$prop)
  se@metadata$TCA_deconv <- res_full
  return(se)
}

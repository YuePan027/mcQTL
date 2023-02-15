#' Estimation of cellular composition in high-throughput data from heterogeneous tissues
#'
#' This function returns a `SummarizedExperiment` object including cell-type proportion estimates for each sample.
#'
#' This is a function developed to implement cell-type proportion deconvolution using either `CIBERSORT` or `nnls`.
#'
#' @param se A `SummarizedExperiment` object with bulk protein/gene expression contained in `counts` slot, and
#' a "signature matrix" which serves as a reference of known cellular signatures contained as an element in `metadata` slot.
#' @param method A character string denotes which deconvolution method to use. In this current version, only `CIBERSORT` or `nnls` is supported.
#'
#' @return A `SummarizedExperiment`. The cell-type proportion estimates for each sample will be stored as an element (`prop`) in `metadata` slot.
#' @export
#'
#' @examples
#'
#' se <- SummarizedExperiment(assays = list(counts = pQTL::protein_data), rowData = pQTL::anno_protein)
#' metadata(se) <- list(sig_matrix = pQTL::ref_data)
#' se <- deconv(se, "cibersort")
#'
#'
deconv <- function(se,
                   method = c("cibersort", "nnls")){
  in_use <- intersect(rownames(assay(se)), rownames(se@metadata$sig_matrix))
  protein_sub <- assay(se)[in_use, , drop=F]
  sig_matrix_sub <- se@metadata$sig_matrix[in_use, , drop=F]

  if(method == "cibersort"){
    result <- CIBERSORT(sig_matrix = sig_matrix_sub,
                        mixture_file = protein_sub,
                        perm=0, QN=TRUE,
                        absolute=FALSE, abs_method='sig.score')
    prop <- result[, 1:ncol(sig_matrix_sub)]
  }

  if(method == "nnls"){
    decon_nnls <- apply(protein_sub, 2, function(y)nnls(sig_matrix_sub,y)$x) %>%t
    prop <- decon_nnls/rowSums(decon_nnls)
    rownames(prop) <- colnames(protein_sub)
    colnames(prop) <- colnames(sig_matrix_sub)
  }

  res_tca <- tca(X = protein_sub,
                 W = prop, refit_W = T,
                 verbose = FALSE)
  prop <- res_tca$W
  se@metadata$prop <- prop
  return(se)
}

#' Estimation of cellular composition in high-throughput data from heterogeneous tissues
#'
#' This function returns a `SummarizedExperiment` object including SNPs used to test for each protein in downstream analysis.
#'
#' This is a function developed to implement cell-type proportion deconvolution using either `CIBERSORT` or `nnls`.
#'
#' @param se A `SummarizedExperiment` object with bulk protein/gene expression contained in `counts` slot.
#' Annotations on each row (protein or gene) should be stored in rowData() with protein or gene symbol as row names
#' The first column should be a character vector indicating which chromosome each protein or gene is on.
#' A "Start" column with numeric values indicating the start position on that chromosome and
#' a "Symbol" column as a unique name for each protein or gene are also required.
#' The information from genetic variants should be stored in a P (the number of SNP) by N (the number of samples, should match the sample in `counts` slot) matrix contained as an element (`SNP_data`) in `metadata` slot.
#' Each matrix entry corresponds to the genotype group indicator (0, 1 or 2) for a sample at a genetic location.
#' The annotations of these SNP should be stored as an element (`anno_SNP`) in `metadata` slot.
#' It should include at least the following columns: "CHROM" (which chromosome the SNP is on),
#' "POS" (position of that SNP) and "ID" (a unique identifier for each SNP, usually a combination of chromosome and its position).
#' @param target_protein A character vector contains proteins/genes names that will be used for downstream analysis.
#' By default, all proteins/genes in `counts` slot will be used.
#' @param target_SNP A character vector contains SNP IDs that will be used for downstream analysis.
#' If not provided, all SNPs will be used for further filtering.
#' @param filter_method A character string denotes which filtering method will be used to filter out unrelated SNPs.
#' If "allele", then the minor allele frequency below argument `filter_allele` will be filtered out.
#' If "distance", then only cis-acting SNPs for each protein (defined as SNPs on the same chromosome and within 1M base pair (bp) range of that protein) will be included for downstream analysis.
#' if "null", then the same SNPs will be used for each protein.
#' @param filter_allele A numeric value denotes the threshold for minor allele frequency.
#'
#' @return A `SummarizedExperiment`. The results after filtering will be stored as an element (`choose_SNP_list`) in `metadata` slot.
#' `choose_SNP_list` is a list with the length of the number of proteins for downstream analysis.
#' Each element stores the index of SNPs to be tested for corresponding protein.
#' The proteins with no SNPs correspond to it will be removed from the returned list.
#'
#' @export
#'
#'
#' @examples
#'
#' se <- SummarizedExperiment(assays = list(counts = pQTL::protein_data),
#'                            rowData = pQTL::anno_protein)
#' metadata(se) <- list(SNP_data = pQTL::SNP_data, anno_SNP = pQTL::anno_SNP)
#' target_protein <- c("NUDT2", "GALT")
#' se <- feature_filter(se, filter_method = c("allele", "distance"), filter_allele = 0.25)
#'
#'
feature_filter <- function(se,
                   target_protein = NULL,
                   target_SNP = NULL,
                   filter_method = c("allele", "distance", "null"),
                   filter_allele = 0.25){

  if (!all(colnames(assay(se)) == colnames(metadata(se)$SNP_data))){
    stop("Samples in protein_data do not match that in SNP_data")
  }

  if (!nrow(metadata(se)$SNP_data) == nrow(metadata(se)$anno_SNP)){
    stop("SNPs contained in annotation data frame `anno_SNP` must match the SNPs in `SNP_data`")
  }

  if(!is.null(target_protein)){
    se <- se[target_protein, ]
  }

  if(!is.null(target_SNP)){
    metadata(se)$SNP_data <- metadata(se)$SNP_data[match(target_SNP, metadata(se)$anno_SNP$ID), , drop=F]
    metadata(se)$anno_SNP <- metadata(se)$anno_SNP[match(target_SNP, metadata(se)$anno_SNP$ID), , drop=F]
  }

  res_TOAST <- NULL
  rowData(se) <- rowData(se)[match(rownames(assay(se)), rowData(se)$Symbol),]
  choose_SNP_list <- rep(list(1:nrow(metadata(se)$SNP_data)), each = nrow(assay(se)))

  if("allele" %in% filter_method){
    prop <- apply(metadata(se)$SNP_data, 1, prop_allele)
    idx <- which(prop > filter_allele & prop <= 1 - filter_allele)
    metadata(se)$SNP_data <- metadata(se)$SNP_data[idx, , drop=F]
    metadata(se)$anno_SNP <- metadata(se)$anno_SNP[idx, , drop=F]
    choose_SNP_list <- rep(list(1:nrow(metadata(se)$SNP_data)), each = nrow(assay(se)))
  }
  if("distance" %in% filter_method){
    choose_SNP_list <- lapply(1:nrow(assay(se)), function(i){
      df <- rowData(se)[which(rowData(se)$Symbol == rownames(assay(se))[i]),]
      idx2 <- which((abs(as.numeric(metadata(se)$anno_SNP$POS) - df$Start) < 1000000) &
                      metadata(se)$anno_SNP$CHROM == df[1,1])
      return(idx2)
    })
    names(choose_SNP_list) <- rownames(assay(se))
    idx3 <- which(unlist(lapply(choose_SNP_list, length)) != 0)
    se <- se[idx3,]
    choose_SNP_list <- choose_SNP_list[idx3]
  }

  metadata(se)$choose_SNP_list <- choose_SNP_list
  return(se)
}

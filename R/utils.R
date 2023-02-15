#' Allele proportion
#'
#' This function calculates the allele frequency.
#'
#' @param x A vector consists of 0, 1 or 2 as genotype group indicator in a genetic location in the genome.
#'
#' @return A numeric value for allele frequency.
#'
#'
prop_allele <- function(x){
  prop <- (sum(x == 0) *2 + sum(x == 1)) / (2*length(x))
  return(prop)
}

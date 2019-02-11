# transform.R

#' MutationTransform
#'
#' This function will transform the data into mutation rate data allowing
#' us to understand the prevalence of mutation by genes and tissue types
#'
#' @param <mut> <Table of all targetted screens by genes, including negatives>.
#' @param <tissue> <boolean, returns mutation rates of genes by tissue if true>.
#' @param <gene> <boolean, returns mutation rates by genes if true>
#' One of gene or tissue must be supplied
#' @return <Matrix of mutation rates>.
#'
#' @author Matthew McNeil


MutationTransform <- function(mut, tissue = T, gene = T) {
  indexList <- list()
  if(tissue){
  indexList$Site <- mut$Site
  }
  if(gene){
  indexList$Gene <- mut$newSymbol
  }
  mutationRates <- tapply(mut$Mutation, indexList, mean)
  return(mutationRates)
}

# [END]




#' Plot network
#'
#' @param go_enrich enrich object
#' @return plot
#' @export
#'
PlotEnrichNetwork <- function(go_enrich=NULL){
    p <- aPEAR::enrichmentNetwork(go_enrich@result, colorBy='pvalue', colorType='pval', nodeSize="Count", drawEllipses=TRUE)

    return(p)
}






#' Calculate gene enrichment by ORA method
#'
#' @param genes gene list
#' @param fgmt gmt file
#' @return enrich object
#' @export
#'
RunORA <- function(
    genes=NULL, 
    fgmt='db.gmt',
    min_n=3,
    max_n=500
){

    db <- clusterProfiler::read.gmt(fgmt)

    enrich_obj <- clusterProfiler::enricher(unique(genes), 
                    TERM2GENE=db, 
                    minGSSize = min_n,
                    maxGSSize = max_n, 
                    pvalueCutoff = 0.1, 
                    qvalueCutoff = 0.2)

    d_ora <- as.data.frame(enrich_obj)

    return(d_ora)
}






#' Calculate gene enrichment by FGSEA method
#'
#' @param data data-frame
#' @param fgmt gmt file
#' @return data frame
#' @export
#'
RunGSEA <- function(
    data=NULL, 
    fgmt='db.gmt',
    min_n=3,
    max_n=500
){
    db <- fgsea::gmtPathways(fgmt)

    data <- data[!duplicated(data$gene),]
    geneList <- data$avg_log2FC
    names(geneList) <- data$gene
    geneList = sort(geneList, decreasing = TRUE)

    d_fgsea <- fgsea::fgseaMultilevel(pathways=db, 
                                stats=geneList,
                                scoreType = "pos",
                                minSize=min_n, 
                                maxSize=max_n)

    d_fgsea <- d_fgsea %>% filter(pval < 0.05) %>% arrange(pval)

    return(d_fgsea)
}
















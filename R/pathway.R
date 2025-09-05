

#' Calculate gene enrichment by ORA method and pathway
#'
#' @param df dataframe
#' @param org_db org db
#'
#' @export
#'
RunPathwayGroup <- function(
    df=NULL, 
    fun='enrichGO',
    species = "Homo sapiens",
    category = "C2",
    org_db=NULL,
    ont='ALL',
    min_n=3,
    max_n=500,
    colors=c('0'='#FFA500', '1'='#4B0082', '2'='#20B2AA'),
    prefix='pathway_group',
    w=8,
    h=5,
    outdir='.'
){
    
    go_enrich <- NULL

    if (fun=='enrichGO'){
        go_enrich <- clusterProfiler::compareCluster(
                        gene~group, 
                        data=df, 
                        fun='enrichGO', 
                        OrgDb=org_db, 
                        keyType = 'SYMBOL', 
                        ont = ont, 
                        readable = T, 
                        pvalueCutoff = 0.01, 
                        minGSSize = 5, 
                        maxGSSize = 500)
    }

    if (fun=='enricher'){
        g_msigdb <- msigdbr::msigdbr(species = species, category = category) %>% dplyr::select(gs_name, gene_symbol)
        go_enrich <- clusterProfiler::compareCluster(
                        gene~group, 
                        data=df, 
                        fun='enricher', 
                        TERM2GENE=g_msigdb,
                        readable = T, 
                        pvalueCutoff = 0.01, 
                        minGSSize = 5, 
                        maxGSSize = 500)
    }

    # First letter to upper case
    firstup <- function(x) {
        substr(x, 1, 1) <- toupper(substr(x, 1, 1))
        x
    }
    go_enrich@compareClusterResult$Description <- firstup(go_enrich@compareClusterResult$Description)
    go_enrich@compareClusterResult$p.adjust <- signif(go_enrich@compareClusterResult$p.adjust, digits=3)
    
    saveRDS(go_enrich, paste0(outdir, "/go_enrich.rds"))

    # 1.output non-redundancy sig. results
    # ONTOLOGY ID Description GeneRatio BgRatio pvalue p.adjust qvalue geneID Count
    write.table(go_enrich@compareClusterResult, file=paste0(outdir, "/ora.", prefix, ".xls"), sep = "\t", row.names = F, quote = F)

    # 2.plot
    d_emap <- enrichplot::pairwise_termsim(go_enrich)
    # 2025-08 version paramater changed
    p <- enrichplot::emapplot(filter(d_emap, pvalue < 0.05 & Count > 1), size_category=.8, pie='Count')
    p <- p + scale_fill_manual(values=colors)


    pdf(paste0(outdir, "/ora.", prefix, ".pdf"), width = 8, height = 5, useDingbats=FALSE)
    print(p)
    dev.off()
}












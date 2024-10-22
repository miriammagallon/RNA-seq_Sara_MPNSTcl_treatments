#########################
#    Packages needed    #
#########################
# library(ReactomePA)
# library(GSEABase)
library(gage)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(magrittr)

# gsea.dir <- paste0("results/GSEA/",celltype)
# genes <- rownames(data.frame(selected_genes$NC$day30))
# sorted.genes <- data.frame(deseq.results$NC$day30[order(deseq.results$NC$day30$log2FoldChange,decreasing = T),])$log2FoldChange
# names(sorted.genes)<-rownames(data.frame(deseq.results$NC$day30[order(deseq.results$NC$day30$log2FoldChange,decreasing = T),]))
# sorted.genes <-na.omit(sorted.genes)
# genes <- length(genes)

#Load annotation data
orgdb <- org.Hs.eg.db
#orgdb <- mouse4302.db
#Enrich functions

#Subset

#x <- enrichGO.BP(gene = genes)

# enrichGO.BP: using a set of genes (ex. DEG), we will obtain enrich BP terms 
enrichGO.BP <- function(genes, orgdb = org.Hs.eg.db) {
  now.msg("Getting enrichemt GO.BP terms...")
  msg("This may take a while...")
  return(enrichGO(gene = genes, OrgDb = orgdb, keyType = "SYMBOL", ont = "BP"))
  
}

# enrichGO.MF: using a set of genes (ex. DEG), we will obtain enrich MF terms 
enrichGO.MF <- function(genes, orgdb = org.Hs.eg.db) {
  now.msg("Getting enrichemt GO.MF terms...")
  msg("This may take a while...")
  return(enrichGO(gene = genes, OrgDb = orgdb, keyType = "SYMBOL", ont = "MF"))  
}

# enrichGO.CC: using a set of genes (ex. DEG), we will obtain enrich CC terms 
enrichGO.CC <- function(genes, orgdb = org.Hs.eg.db) {
  now.msg("Getting enrichemt GO.CC terms...")
  msg("This may take a while...")
  return(enrichGO(gene = genes, OrgDb = orgdb, keyType = "SYMBOL", ont = "CC"))  
}


# enrichReactome.Paths: using a set of genes (ex. DEG), we will obtain enrich Reactome pathways terms 
enrichReactome.Paths <- function(genes, orgdb = org.Hs.eg.db){
  now.msg("Getting enrichemt Reactome pathways terms...")
  msg("This may take a while...")
  #We need to change gene symbols to Entrez ID.
  symb.to.entrezid <- bitr(geneID = genes, OrgDb = orgdb ,fromType='SYMBOL', toType='ENTREZID',drop = TRUE)
  symb.to.entrezid[is.na(symb.to.entrezid)] <- NA
  genes <- symb.to.entrezid$ENTREZID
  
  enrich.react.paths <- clusterProfiler::enrichKEGG(gene = genes,
                                                    organism = "hsa",
                                                    keyType = "kegg",
                                                    minGSSize = 15)
  enrich.react.paths <- DOSE::setReadable(x = enrich.react.paths, OrgDb = orgdb, keytype = "ENTREZID")
  
  now.msg("Enrichment Reactome Pathway terms obtained")
  return(enrich.react.paths)
}


# enrichKEGG.Paths: using a set of genes (ex. DEG), we will obtain enrich KEGG pathways terms 
enrichKEGG.Paths <- function(genes, orgdb = org.Hs.eg.db ){
  now.msg("Getting enrichemt KEGG pathways terms...")
  msg("This may take a while...")
  #We need to change gene symbols to Entrez ID.
  symb.to.entrezid <- bitr(geneID = genes, OrgDb = orgdb ,fromType = 'SYMBOL', toType = 'ENTREZID',drop = TRUE)
  symb.to.entrezid[is.na(symb.to.entrezid)] <- NA
  genes <- symb.to.entrezid$ENTREZID
  
  enrich.KEGG <- clusterProfiler::enrichKEGG(gene = genes,
                                             organism = "hsa",
                                             keyType = "kegg",
                                             minGSSize = 15)
  # enrich.KEGG <- DOSE::setReadable(x = enrich.KEGG, OrgDb = orgdb, keyType =  "ENTREZID")
  enrich.KEGG <- DOSE::setReadable(x = enrich.KEGG, OrgDb = orgdb, keyType =  "ENTREZID")
  now.msg("Enrichment KEGG Pathway terms obtained")
  return(enrich.KEGG)
}


#All enrichmed term list
getAllEnrichTermList <- function(genes){
  now.msg("Getting all enriched terms...")
  msg("This may take a while...")
  enrichtermtlist <- list(BP = enrichGO.BP(genes),
                          MF = enrichGO.MF(genes),
                          CC = enrichGO.CC(genes))#),
                        #  KEGG = enrichKEGG.Paths(genes))
  
  
  now.msg("All enriched terms obtained")
  return(enrichtermtlist)
  
}


##################
#  Sorted genes  #
##################

# getGSEAGO.BP: using all genes to obtain gene enrichment by sorting all genes by LFC, GO-BP. 
getGSEAGO.BP <- function(sorted.genes, qvalue.cutoff = 0.05){
  now.msg("Running GSEA GO.BP...")
  msg("This may take a while...")
  return(clusterProfiler::gseGO(sorted.genes,
                                ont = "BP",
                                OrgDb = orgdb,
                                keyType = "SYMBOL",
                                pvalueCutoff = qvalue.cutoff,
                                minGSSize = 15,
                                pAdjustMethod = "fdr",
                                by = "fgsea"))
  now.msg("GSEA GO.BP done")
  
}


# getGSEAGO.MF: using all genes to obtain gene enrichment by sorting all genes by LFC, GO-MF. 
getGSEAGO.MF <- function(sorted.genes, qvalue.cutoff = 0.05){
  now.msg("Running GSEA GO.MF...")
  msg("This may take a while...")
  return(clusterProfiler::gseGO(sorted.genes,
                                ont = "MF",
                                OrgDb = orgdb,
                                keyType = "SYMBOL",
                                pvalueCutoff = qvalue.cutoff,
                                minGSSize = 15,
                                pAdjustMethod = "fdr",
                                by = "fgsea"))
  now.msg("GSEA GO.MF done")
}


# getGSEAGO.CC: using all genes to obtain gene enrichment by sorting all genes by LFC, GO-CC. 
getGSEAGO.CC <- function(sorted.genes, qvalue.cutoff = 0.05){
  
  now.msg("Running GSEA GO.CC...")
  msg("This may take a while...")
  return(clusterProfiler::gseGO(sorted.genes,
                                ont = "CC",
                                OrgDb = orgdb,
                                keyType = "SYMBOL",
                                pvalueCutoff = qvalue.cutoff,
                                minGSSize = 15,
                                pAdjustMethod = "fdr",
                                by = "fgsea"))
  now.msg("GSEA GO.CC done")
  
}
# 
# xx <- getGSEAGO.BP(sorted.genes = sorted.genes)
# head(data.frame(xx))


# getGSEAReactomePathways: using all genes to obtain gene enrichment by sorting all genes by LFC, ReactomePathways. 
getGSEAReactomePathways <- function(sorted.genes, qvalue.cutoff = 0.05){
  now.msg("Running GSEA Reactome Pathways...")
  msg("This may take a while...")
  symb.to.entrezid <- bitr(geneID = names(sorted.genes),
                           OrgDb = orgdb,
                           fromType='SYMBOL',
                           toType='ENTREZID',
                           drop = TRUE)
  
  symb.to.entrezid <- symb.to.entrezid[!duplicated(symb.to.entrezid$SYMBOL),]
  names(sorted.genes)<- symb.to.entrezid$ENTREZID
  sorted.genes <- sorted.genes[!is.na(names(sorted.genes))]
  
  gsea.reactome.paths <- ReactomePA::gsePathway(geneList = sorted.genes, minGSSize = 15)
  gsea.reactome.paths <- DOSE::setReadable(gsea.reactome.paths, OrgDb = orgdb, keyType ="ENTREZID")
  
  return(gsea.reactome.paths)
  now.msg("GSEA GSEA Reactome Pathways done")
}


# getGSEAKEGGpathways: using all genes to obtain gene enrichment by sorting all genes by LFC, KEGG Pathways. 
getGSEAKEGGpathways <- function(sorted.genes, qvalue.cutoff = 0.05){
  now.msg("Running GSEA KEGG Pathways...")
  msg("This may take a while...")
  sorted.genes <- sorted.genes[!is.na(names(sorted.genes))]
  symb.to.entrezid <- bitr(geneID = names(sorted.genes),
                           OrgDb = orgdb,
                           fromType='SYMBOL',
                           toType='ENTREZID',
                           drop = TRUE)
  
  symb.to.entrezid <- symb.to.entrezid[!duplicated(symb.to.entrezid$SYMBOL),]
  
  names(sorted.genes)<- symb.to.entrezid$ENTREZID
  sorted.genes <- sorted.genes[!is.na((sorted.genes))]
  
  gsea.KEGG.pathways <- clusterProfiler::gseKEGG(geneList = sorted.genes,
                                                 organism = "hsa",
                                                 keyType = "kegg",
                                                 minGSSize = 15 )
  
  gsea.KEGG.pathways <- DOSE::setReadable(gsea.KEGG.pathways,
                                          OrgDb = orgdb,
                                          keyType ="ENTREZID")
  
  return(gsea.KEGG.pathways)
  now.msg("GSEA GSEA KEGG Pathways done")
}

fromSymbolToEntrezID <- function(sorted.genes){
  sorted.genes <- sorted.genes[!is.na(names(sorted.genes))]
  symb.to.entrezid <- bitr(geneID = names(sorted.genes),
                           OrgDb = orgdb,
                           fromType='SYMBOL',
                           toType='ENTREZID',
                           drop = TRUE)
  
  symb.to.entrezid <- symb.to.entrezid[!duplicated(symb.to.entrezid$SYMBOL),]
  names(sorted.genes)<- symb.to.entrezid$ENTREZID
  sorted.genes <- sorted.genes[!is.na((sorted.genes))]
  return(sorted.genes)
}

# kegg <- getGSEAKEGGpathways(sorted.genes)
# 
# enrich.or.GSEA.object <- kegg
# enrich.or.GSEA.object <- x

# selectSpecificTerms: Once we he all enriched terms or GSEA output, we might decide to plot only specific genes.
# To get it, we have developed one function --> selectSpecificTerms.

selectSpecificTerms <- function(enrich.or.GSEA.object, GO.KEGG.reactome.term){
  
  selected.terms <- na.omit(enrich.or.GSEA.object[GO.KEGG.reactome.term])
 
  # We must build an S4 object in order to be able to represent specific terms in cnetplot or emaplot.
  if(class(enrich.or.GSEA.object) == "gseaResult"){
    setClass("gseaResult", representation = representation(result = "data.frame", 
                                                           params = "list",
                                                           organism = "character",
                                                           setType = "character",
                                                           keytype = "character", 
                                                           readable = "logical",
                                                           geneList = "numeric",
                                                           gene2Symbol = "character"))
    selected.terms <- new("gseaResult", result = selected.terms)
    selected.terms@params <- enrich.or.GSEA.object@params
    selected.terms@organism <- enrich.or.GSEA.object@organism
    selected.terms@setType <- enrich.or.GSEA.object@setType
    selected.terms@keytype <- enrich.or.GSEA.object@keytype
    selected.terms@readable <- enrich.or.GSEA.object@readable
    selected.terms@geneList <- enrich.or.GSEA.object@geneList
    selected.terms@gene2Symbol <- enrich.or.GSEA.object@gene2Symbol
    
  }else if (class(enrich.or.GSEA.object) =="enrichResult"){
    
    
    # We must build an S4 object in order to be able to represen specific terms in cnetplot or emtplot.
    setClass("enrichResult", representation = representation(result = "data.frame",
                                                             pvalueCutoff = "numeric",
                                                             qvalueCutoff = "numeric",
                                                             pAdjustMethod = "character",
                                                             organism = "character",
                                                             ontology = "character",
                                                             readable = "logical",
                                                             gene = "character",
                                                             gene2Symbol = "character",
                                                             keytype = "character"))
    
    selected.terms <- new("enrichResult", result = selected.terms,
                          pvalueCutoff = enrich.or.GSEA.object@pvalueCutoff,
                          qvalueCutoff = enrich.or.GSEA.object@qvalueCutoff,
                          pAdjustMethod = enrich.or.GSEA.object@pAdjustMethod,
                          organism = enrich.or.GSEA.object@organism,
                          ontology = enrich.or.GSEA.object@ontology,
                          readable = enrich.or.GSEA.object@readable,
                          gene = enrich.or.GSEA.object@gene,
                          gene2Symbol = enrich.or.GSEA.object@gene2Symbol,
                          keytype = enrich.or.GSEA.object@keytype) 
    
  }
  
  return(selected.terms)
}


#################
#     Plots     #
#################

# x
# specific.terms <- c("GO:0045047","GO:0006613")
# enrich.GSEA.terms <- x
# head(data.frame(kegg),5)
# s <- selectSpecificTerms(enrich.or.GSEA.object = kegg, GO.KEGG.reactome.term =  specific.terms)

# plotCNET: plot genes associated with its terms as well as its expresion in GSEA terms.
plotCNET <- function(enrich.GSEA.terms, specific.terms, sorted.genes = NULL, select.specific.terms = FALSE, num.terms.to.show = 5, title, title.size =20, lab.size=20){
  
  if(select.specific.terms == TRUE){
    #If select.specific.terms is TRUE we will select those terms specified in specific term parameters before represent.
    selected.terms <- selectSpecificTerms(enrich.GSEA.terms, GO.KEGG.reactome.term = specific.terms)
    
    cnet <- cnetplot(selected.terms, showCategory = nrow(selected.terms), foldChange = sorted.genes) + 
      ggtitle(label = title)+ theme(text = element_text(size=lab.size), title =element_text(size = title.size))
    print(cnet)
    
    
  }else{
    cnet <- cnetplot(enrich.GSEA.terms, showCategory = num.terms.to.show, foldChange = sorted.genes ) + 
      ggtitle(label = title)+ theme(text = element_text(size=lab.size), title =element_text(size = title.size))
   print(cnet)
  }
  
  invisible(enrich.GSEA.terms)
}
# sorted.genes <- 0
# specific.terms <- s
# plotEMAPPLOT(enrich.GSEA.terms = cc)


# plotEMAPPLOT: This function will plot the association between terms. We can also select specific terms and we can observe,
# its assocation.
plotEMAPPLOT <- function(enrich.GSEA.terms, select.specific.terms = FALSE, specific.terms, title, title.size =20, lab.size = 20){
  
  if(select.specific.terms == TRUE){
    #If select.specific.terms is TRUE we will select those terms specified in specific term parameters before represent.
    selected.terms <- selectSpecificTerms(enrich.GSEA.terms, GO.KEGG.reactome.term = specific.terms)
    print(emapplot(selected.terms, showCategory = nrow(selected.terms)) + ggtitle(label = title)+
            theme( title = element_text(size=title.size), text = element_text(size=lab.size)))
    
    
  }else{
   print(emapplot(enrich.GSEA.terms) + 
           ggtitle(label = title)+
           theme( title = element_text(size=title.size), text = element_text(size=lab.size)))
  }
  invisible(enrich.GSEA.terms)
}


# plotHeatplot: this function is used to represent how ranked genes are associeated between different terms.
plotHeatplot <- function(GSEA.terms, sorted.genes, num.terms.to.show = 5, title, title.size = 20, lab.size = 20){
  print(heatplot(x = GSEA.terms, showCategory = num.terms.to.show, foldChange = sorted.genes)
          + ggtitle(label = title) + theme(text = element_text(size =lab.size),axis.title = element_text(size =title.size)) 
          )
  invisible(GSEA.terms)
}


# plotDotplot: this function is used to represent enriched terms or GSEA terms as dots.
plotDotplot <- function(enrich.GSEA.terms, title, lab.size = 20, tilte.size = 10){
  print(dotplot(enrich.GSEA.terms)+ ggtitle(label = title) + 
          theme(text = element_text(size =lab.size),axis.title = element_text(size =title.size))) 
  invisible(enrich.GSEA.terms)
}


# plotGSEAplot: this function is used to represent GSEA terms as Broadig institute GSEA classic plots.
plotGSEAplot <- function(GSEA.terms, geneSetID){
  print(gseaplot(GSEA.terms, geneSetID = geneSetID, title = GSEA.terms@result[geneSetID,"Description"]))
  invisible(GSEA.terms)
}


# We can also represent GSEA enriched pathways using pathview funciton.
# pathview(gene.data = sorted.genes, pathway.id = "hsa05230")
# 
# 
# selectSpecificTerms(x, GO.KEGG.reactome.term = GO.KEGG.reactome.term)%>%plotCNET()%>%plotDotplot()
# plotCNET(x)
# 
# enrichGO.BP(genes) %>% plotCNET(select.specific.terms = FALSE)
# bp <-enrichGO.BP(names(lfc.gene))
# head(data.frame(bp),25)
# kegg <- enrichKEGG.Paths(names(lfc.gene))
# plotCNET(bp)
# gdf <- deseq.results$NC$day30
# lfc <-gdf[order(gdf$log2FoldChange,decreasing = T),]
# lfc.gene <- lfc$log2FoldChange
# names(lfc.gene) <- rownames(lfc)
# entrezid <- bitr(geneID = names(lfc.gene),
#                  OrgDb = org.Hs.eg.db,
#                  fromType='SYMBOL',
#                  toType='ENTREZID')
# entrezid <- entrezid[which(!duplicated(entrezid$SYMBOL)),]
# # 
# nrow(entrezid)
# # 
# # 
# names(lfc.gene) <- (entrezid$ENTREZID)
# paths <- c("hsa04360", "hsa04390")
# # 
# for (i in seq_len(length(paths))){
#   pathview(lfc.gene, pathway.id = paths[i], kegg.dir = "results/pathways",out.suffix = "NCvsday30")
# }

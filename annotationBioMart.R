
library(biomaRt)

# ensembl=useMart("ensembl")
# listDatasets(ensembl)
# getBiomartGOAnnotation(genes = k)
# ensembl=useMart("ensembl")
# ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
# listMarts(ensembl)
#' getBioMartGOAnnotation
#'
#' @description
#'
#' getBioMartGOAnnotation function provides  GO terms annatotion of specific genes.
#'
#' @details
#'
#' getBioMartGOAnnotation function provides  GO terms annatotion of specific genes using BioMart annotation.
#' The user simply enters the genes that want to annote, what gene identification
#' the samples have, and the reference genome. This information is known using listFilter() and
#' listDataSet(), respectively.
#'
#' @usage
#'
#' getBioMartGOAnnotation (genes, gene_id = 'hgnc_symbol',
#' ref_genome = "hsapiens_gene_ensembl" ,  verbose = TRUE )
#'
#' @param genes (character) character vector containing all genes that wanted to be analysed.
#' @param gene_id (character) gene identifications. To know which gene identification has to be used,
#' the user can use listFilter() of biomaRt package.(default: gene_id = 'hgnc_symbol')
#' @param ref_genome (character) reference genome of the sample. To know which genome to use,
#' the user can us listDataSet() of biomaRt package. (default: ref_genome = "hsapiens_gene_ensembl")
#' @param verbose (logical) whether the user wants some messages over the process.(default: verbose = TRUE)
#'
#' @note Sometimes BioMart web service is down because of maintenance. It migth be this function does not
#' work for that reason.
#'
#' @return
#'
#' getBioMartGOAnnotation's main output is a list of data frames classifaying the GO terms annotation code,
#' and information about the GO term annoted by the biological function they have (CC, MF or BP).
#'
#' @export getBioMartGOAnnotation

getBioMartGOAnnotation <- function(genes, gene_id = 'hgnc_symbol', 
                                   ref_genome = "hsapiens_gene_ensembl" ,  verbose = TRUE ){
  now.msg("Getting Biomart GO annotation...",verbose = verbose)
  ensembl <- useMart(biomart = "ensembl", dataset = ref_genome)
  
  annot <- getBM(attributes=c("hgnc_symbol", 'entrezgene_id', 'go_id', "name_1006", "namespace_1003"),
      filters = gene_id, 
      values = genes, 
      mart = ensembl)
    
  BP <- data.frame(annot)[which(annot$namespace_1003 == "biological_process"),]
  MF <- data.frame(annot)[which(annot$namespace_1003 == "molecular_function"),]
  CC <- data.frame(annot)[which(annot$namespace_1003 == "cellular_component"),]
  data <- list(BP =data.frame(BP), MF = data.frame(MF), CC = data.frame(CC))
  
  now.msg("Biomart GO annotation done",verbose = verbose)
  return(data)
}



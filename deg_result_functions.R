#####################################################################################
#                    Differencial Expressed Genes Functions                         #  
#####################################################################################

source(file ="./RNA-seq_Sara_MPNSTcl_treatments/utils.R" )

#sample_group == column from samples data frame

#' getFilteredDDS
#'
#' @description
#'
#' getFilteredDDS function gets a filtered DESeqDataSet object where the design formula is added depending on
#' the test performed.
#'
#' @details
#'
#' getFilteredDDS function gets a filtered DESeqDataSet object. This funtion contains DESeqDataSetFromTximport function
#' from DESeq2 package, used to store the input values from tximport list. Depending on the test the user wants to do,
#' the design formula of DESeqDataSetFromTximport will vary. If the user like studying more than one sample group, LRT test
#' will be performed. If the interaction between different sample groups is analysed, then this function will join in
#' a unique sample group variable as one condition to study. If the user like to study one condition (sample_group parameter), it will
#' be obtained a simple formula design using any test.
#'
#' @usage
#'
#' getFilteredDDS(tximport, samples_group, samples_df, filter_min_reads = 5,
#' filter_min_samples  = 1, test = "Wald", interaction = FALSE, verbose = TRUE)
#'
#' @param tximport (list) tximport object produced by tximport function
#' @param samples_group (character) name of the factor from samples_df to contrast their levels. If sample_group is
#' more than one, LRT  test must be performed.
#' @param samples_df (data.frame) data frame containing all information about the samples.
#' @param filter_min_reads (numeric) minimum of count reads to filter (default: filter_min_reads = 5)
#' @param filter_min_samples (numeric)minimun of columns to apply filter_min_reads(default: filter_min_samples = 1).
#' @param test (character) test that is going to be performed in order to get DESeq design formula.
#' it could be whether "Wald" or "LRT" tests (default: test = "Wald").
#' @param interaction (logical) only used if LRT test is performed and sample_group are two.
#' Whether the design formula study the interaction of sample_group. (default: ineraction = FALSE)
#' @param verbose (logical) whether the user like some messages over the process.
#'
#'
#' @return
#'
#' getFilteredDDS's main output is a filtered DESeqDataSet object ready to get DEG.
#'
#' @export getFilteredDDS
#'
#'
getFilteredDDS <- function(tximport, samples_group, samples_df, filter_min_reads = 5,
                           filter_min_samples  = 1, test = "Wald", interaction = FALSE,
                           verbose = TRUE){
  
  now.msg(" Filtering DESeqDataSet...", verbose = verbose)
  
  if(test == "Wald"){
    
    if(length(samples_group) == 1){
      
      # Generation of a factor to contrast
      cont<-samples_df[,samples_group]
      df<- cbind(samples_df, design = as.character(cont), stringsAsFactors=FALSE)
      
      #DESFromTXI
      dds <- DESeqDataSetFromTximport(txi = tximport, colData = df, design = ~ design)
      
      
    }else{
      stop("Wald test is only for one conditon and sample_group's parameter has more than one")
    }
  }
  
  if(test == "LRT"){
    
    msg(" generating LRT test design for contrast  ...", verbose = verbose)
    
    # Generation of a factor to contrast
    
    if (samples_group > 1 && interaction == FALSE){
     
       # Generation of a factor to contrast
      
      cont<-samples_df[,samples_group]
      df<- cbind(samples_df, design = cont, stringsAsFactors=FALSE)
      
      #DESFromTXI
      
      for( i in 1:(length(samples_group)-1)){
        dds <- DESeqDataSetFromTximport(txi = tximport, colData = df,
                                        design =formula(paste0("~","design.",samples_group[i],
                                                               "+", "design.", samples_group[i+1])))
      }
      
    }else if ( samples_group == 1 && interaction == FALSE){
      # Generation of a factor to contrast
      
      cont<-samples_df[,samples_group]
      df<- cbind(samples_df, design = as.character(cont), stringsAsFactors=FALSE)
        
      #DESFromTXI
      dds <- DESeqDataSetFromTximport(txi = tximport, colData = df, design = ~ design)
        
      
      
    }else if (samples_group == 1 && interaction == TRUE){
      
      msg(" Preparing LRT design for interaction contrast...", verbose = verbose)
      
      for( i in 1:(length(samples_group)-1)){
      df$design = factor(paste0(df[,paste0("design.",samples_group[i])],"_",
                                df[,paste0("design.",samples_group[i+1])]))
      }
      
      dds <- DESeqDataSetFromTximport(txi = tximport, colData = df, 
                                      design = ~design)
    }
  }
  
  #estimateSizeFactor
  dds <- estimateSizeFactors(dds)
  
  # Filtering data
  msg("Fitering data...", verbose = verbose)
  keep <- rowSums(counts(dds) >= filter_min_reads) >= filter_min_samples
  filtered_dds <- dds[keep,]
  
  now.msg(" DESeqDataSet Filtered", verbose = verbose)
  
  return(filtered_dds)
  
}

#' runDESeq
#'
#' @description
#'
#' runDESeq function gives a DESeqResults object containing all DEG ordered by what the user decides ("padj" value
#' or "log2FoldChange") only specifying the name of the samples that want to be analysed, the reference group of the contrast,
#' the DESeqDataSet object filtered previously, the data frame where sample information is contained,
#' and the name of the variable in the sample information where sample names are annotated.
#'
#'
#' @details
#'
#' runDESeq function gets a filtered DESeqDataSet object. This function can get DEG only specifying the name of the samples that want to be analysed, the reference group of the contrast,
#' the DESeqDataSet object filtered previously, the data frame where sample information is contained,
#' and the name of the variable in the sample information where sample names are annotated.
#' This funtion contains DESeqDataSetFromTximport function from DESeq2 package, used to store the input values
#' from tximport list. Depending on the test that the user wants to do ("LRT" or "Wald"),
#' the design formula of DESeqDataSetFromTximport() will be modify. If the user analyses more than one sample group,
#' LRT test must be performed. If the interaction between different sample groups is analysed (interaction = TRUE),
#' then, this function will make the reduce model of de disign acquired with getFilteredDDS()function.
#' (NOTE: if design formula wants to be known, use design() formula over DESeqDataSet).
#' If two samples_group  are analysed without interaction (interaction = FALSE), the user must use studied_condition
#' parameter and give the name of the sample_group (condition) from which  the user wants to obtain the DEG.
#' By this way, runDESeq function will build the reduce design formula and will get de DEG which will show any difference
#' across the sample_group which was specified in studied_condition.
#' If the user only wants to get DEG of one condition, the test used might be "Wald" or "LRT". This function will construct
#' the reduce model of the design formula if LRT is performed, and will get DEG ordered by ordered_by parameter (padj or log2FoldChange).
#' The control group of the sample_group levels of any type of contrast must be specified in ref_group parameter.
#' Moreover, this function has the possibility to perform two different lfcShrink test to get de DEG; the "normal" one
#' which is the default in this function, and "apeglm". If the number of the samples and replicates are small,
#' "normal" lfcShrink analysis is recomended.
#'
#' @usage
#'
#' runDESeq (sample_names, ref_group, filtered_dds, samples_df, sname_variable,
#' type_lfcShrink = "normal", maxpvalue = 0.05, maxlfc = 0, order_by = "padj",
#' test = "Wald", studied_condition = NULL, interaction = FALSE,  verbose = TRUE, ...)
#'
#' @param sample_names (character) name of the samples that the user wants to analyse.
#' @param ref_group (character) the level control group of the contrast.
#' @param filtered_dds (DESeqDataSet) the DESeqDataSet object obtained by getFilteredDDS()function.
#' @param samples_df(data.frame) data frame containing all information about the samples
#' @param sname_variable (character) name of the variable where sample_names is contained in samples_df.
#' @param type_lfcShrink (character) whether to do "normal" lfcShrik or "apeglm" lfcShrink (default:type_lfcShrink = "normal").
#' @param maxpvalue (numeric) maximum p-value to get DEG (default: maxpvalue = 0.05)
#' @param maxlfc (numeric) maximum log2FoldChange to get DEG (default: maxlfc = 0)
#' @param order_by (character) whether to order the results by "padj" or "log2FoldChange" (default: order_by = "padj")
#' @param test (character) whether to do "Wald" test or "LRT" (default:test = "Wald").
#' @param studied_condition (character) name of the sample_group
#' from which the user wants to obtain DEG, when interaction = FALSE. (default:studied_condition = NULL)
#' @param interaction (logical) whether the interaction between sample_groups (more than one) is studied. NOTE: LRT test must be performed.
#' (default: interaction = FALSE)
#' @param verbose(logical) whether the user like some messages over the process (default: verbose = TRUE).
#' @param ... further arguments passed to getFilteredDDS
#'
#'
#' @return
#'
#' runDESeq's main output is a DESeqResults object containing DEG ordered by padj or log2FoldChange.
#'
#' @export runDESeq
#'

runDESeq <- function(sample_names, ref_group, filtered_dds, samples_df, 
                      sname_variable, type_lfcShrink = "normal", maxpvalue = 0.05, 
                      maxlfc = 0, order_by = "padj", test = "Wald", studied_condition = NULL, 
                      interaction = FALSE,  verbose = TRUE, ...) {
  
  now.msg(" Starting DESeq...", verbose = verbose)
  
  # Generation of a factor to contrast
  
  cont<-data.frame(colData(filtered_dds)[, grepl( "design",names(colData(filtered_dds)))])
  sample_group_levels <- cont[, grepl( "design",names(cont))][samples_df[,sname_variable] %in% sample_names]
  sample_group_levels <- unique(sample_group_levels) 
  
  
  #DESeq
  
  if(test == "LRT"){
    
    msg(" Calculating Likelyhood Ratio Test DESeq ", verbose = verbose)
    
    if(interaction == FALSE && length(as.data.frame(sample_group_levels)) == 1){
      column <- "design"
      colData(filtered_dds)[,column] <- droplevels(relevel(colData(filtered_dds)[,column], 
                                                           ref = ref_group))
      dd <- DESeq(filtered_dds, full = design(filtered_dds), test = "LRT", reduced = ~1, ...)
      other_group <- as.character(sample_group_levels[!(sample_group_levels %in% ref_group)])
      
      
    }else if(interaction == FALSE && length(as.data.frame(sample_group_levels)) > 1){
      
      column <- paste0(colnames(cont)[colnames(cont)!=paste0("design.",studied_condition)])
      colData(filtered_dds)[,column] <- droplevels(relevel(colData(filtered_dds)[,column], 
                                                           ref = ref_group))
      dd <- DESeq(filtered_dds, full = design(filtered_dds), test = "LRT",
                  reduced = formula(paste0("~", column)), ...)
      other_group <- as.character(sample_group_levels[,column][!(sample_group_levels[,column] %in% ref_group)])
      
      
    }else if(interaction == TRUE && length(as.data.frame(sample_group_levels)) > 1){
      column <- "design"
      ref_group <- as.character(sample_group_levels[,column][grepl( ref_group,sample_group_levels[,column])])
      colData(filtered_dds)[,column] <- droplevels(relevel(colData(filtered_dds)[,column], 
                                                           ref = ref_group))
      dd <- DESeq(filtered_dds, full = design(filtered_dds), test = "LRT", reduced = ~1, ...)
      other_group <- as.character(sample_group_levels[,column][!(sample_group_levels[,column] %in% ref_group)])
    }
    
  }else if( test == "Wald"){
    
    msg(" Calculating Wald Test DESeq ", verbose = verbose)
    column <- "design"
    colData(filtered_dds)[,column]<- droplevels(relevel(colData(filtered_dds)[,column], ref = ref_group))
    dd <- DESeq(filtered_dds, full = design(filtered_dds), test = "Wald", ...)
    other_group <- as.character(sample_group_levels[!(sample_group_levels %in% ref_group)])
    #dd<- DESeq(filtered_dds, full = design(filtered_dds), test = "Wald")
    
  }else{
    stop("Test parameter must be Wald or LRT")
  }
  
  msg(" Getting DESeq results ", verbose = verbose)
  #other_group <- as.character(sample_group_levels[!(sample_group_levels %in% ref_group)])
  
  if (type_lfcShrink == "apeglm"){
    
    msg("Approximate Posterior Estimation for  the  GLM (apeglm) lfcShrink ")
    
    if(test == "Wald"){
      
      res <- results(dd, contrast = c( column, other_group, ref_group ), alpha = maxpvalue,
                     lfcThreshold = maxlfc)
      res <- lfcShrink(dd, res = res, coef = paste0(column,"_",other_group, "_vs_", ref_group), 
                       type= type_lfcShrink)
      res <- res[order(res[,order_by]),]
      
    }else{
      
      res <- results(dd, contrast = c( column, other_group, ref_group ), alpha = maxpvalue)
      res <- lfcShrink(dd, res = res, coef = paste0(column,"_",other_group, "_vs_", ref_group),
                       type= type_lfcShrink)
      res <- res[order(res[,order_by]),]
      
    }
    
  } else if(type_lfcShrink == "normal"){
    msg("Normaly-distributed prior (Normal) lfcShrink ")
    
    if(test == "Wald"){
      
      res <- results(dd, contrast = c( column, other_group, ref_group ), alpha = maxpvalue, 
                     lfcThreshold = maxlfc)
      res <- lfcShrink(dd, res = res, contrast = c( column, other_group, ref_group ), 
                       type= type_lfcShrink)
      res <- res[order(res[,order_by]),]
      
    }else{
      
      res <- results(dd, contrast = c( column, other_group, ref_group ), alpha = maxpvalue)
      res <- lfcShrink(dd, res = res, contrast = c( column, other_group, ref_group ), 
                       type= type_lfcShrink)
      res <- res[order(res[,order_by]),]
      
    }
  }
  
  
  now.msg("DESeq done", verbose = verbose)
  return (res)
}


#' getTopGenesLFC
#'
#' @description
#'
#' getTopGenesLFC function gets from a DESeqResults object only those genes with a padj_value gave by the user,
#' and ordered by LFC.
#'
#'
#' @details
#'
#' getTopGenesLFC function gets from a DESeqResults objet only those genes with a padj_value gived by the user,
#' and ordered by absolute LFC.
#'
#' @usage
#'
#' getTopGenesLFC(deseq_results, padj_value)
#'
#' @param deseq_results (character) name of the samples that the user wants to analyse.
#' @param padj_value (numeric) padjusted vaule to select genes (default: padj_value = 0.05)
#'
#' @return
#'
#' getTopGenesLFC's main output is a data frame only with those genes selected by padj_value and ordered by LFC.
#'
#' @export getTopGenesLFC
#'

getTopGenesLFC <- function (deseq_results, padj_value = 0.05){
  
  deseq_results<- data.frame(deseq_results)
  
  # We select those genes which have a pvalue <= padj_value from all analised genes.
  # Once we have selected those genes, we proceed to order them by log2FoldChange absolute value (decreasing).
  # By this way, we will have the most upregulated or downregulated genes at the begining of the table.
  
  selected_genes<- deseq_results[deseq_results$padj<= padj_value,]
  selected_genes<- selected_genes[order(abs(selected_genes$log2FoldChange),decreasing = T),][which(selected_genes$log2FoldChange != "NA"),]
  return(selected_genes)
}

#################
getGroupSpecificGenes <- function(select_genes, query_group, updown_markers){
  #Preparing the data to analyse
  
  #ref_list_group will have information of query_group as contrast reference
  
  ref_list_group <- select_genes[[query_group]]
  names(ref_list_group)
  #query_list_group will have information of query_group as contrast query
  
  query_list_group<-list()
  for( i in 1:length(select_genes)){
    ref <- names(select_genes[i])
    query_list_group[[ref]] <- select_genes[[i]][[query_group]]
  }
  
  if( updown_markers == "up"){
    #By this way over expressed genes from query_group are obtained
    
    select_ref_up <- list()
    if(length(ref_list_group )>0){
      for(i in 1:length(ref_list_group)){
        nm<- names(ref_list_group[i])
        select_ref_up[[nm]] <- ref_list_group[[nm]][(ref_list_group[[nm]]["log2FoldChange"])< 0,][which(ref_list_group[[nm]]["log2FoldChange"]!="NA") ,]
      }
    }else{
      select_ref_up <- select_genes[[query_group]]
    }
    
    names(select_ref_up)
    selectq_up<- list()
    if(length(query_list_group)>0){
      
      for(i in 1:length(query_list_group)){
        nmq<- names(query_list_group[i])
        selectq_up[[nmq]]<- query_list_group[[nmq]][(query_list_group[[nmq]]["log2FoldChange"])>0,][which(query_list_group[[nmq]]["log2FoldChange"]!="NA") ,]
      }
      
    }else{
      selectq_up<-query_list_group
      
    }
    #Selection of common genes between all contrast made
    data_up<- c(selectq_up, select_ref_up)
    common_names = Reduce(intersect, lapply(data_up, row.names))
    data_up = lapply(data_up, function(x) { x[row.names(x) %in% common_names,] })
    data_up <- lapply(data_up, function(x) x[!is.na(x$log2FoldChange),])
    
    #we order stage specifc markers regarding LFC mean of the gene in the different stages in a decrasing way
    data <- seq_len(nrow(data.frame(data_up[1])))
    names(data) <- rownames(data.frame(data_up[1]))
    
    for (g in names(data)){
      data[g] <- mean(unlist(lapply(data_up,function(x) abs(x[g,"log2FoldChange"]))))
      
    }
    data <- sort(data, decreasing = T)
    data <- names(data)
    
    
    
    
  }else if(updown_markers == "down"){
    #By this way infra-expressed genes from query_group are obtained
    
    select_ref_down <- list()
    if(length(ref_list_group )>0){
      
      for(i in 1:length(ref_list_group)){
        nm<- names(ref_list_group[i])
        select_ref_down[[nm]] <- ref_list_group[[nm]][(ref_list_group[[nm]]["log2FoldChange"])>0,]#[which(ref_list_group[[nm]]["log2FoldChange"]!="NA") ,]
      }
      
    }else{
      select_ref_down <- ref_list_group
    }
    
    selectq_down<- list()
    if(length(query_list_group )>0){
      
      for(i in 1:length(query_list_group)){
        nmq<- names(query_list_group[i])
        selectq_down[[nmq]]<- query_list_group[[nmq]][(query_list_group[[nmq]]["log2FoldChange"])<0,][which(query_list_group[[nmq]]["log2FoldChange"]!="NA") ,]
      }
      
    }else {
      selectq_down <- query_list_group
      
    }
    #Selection of common genes between all contrast made
    data_down<- c(selectq_down, select_ref_down)
    common_names = Reduce(intersect, lapply(data_down, row.names))
    data_down <-lapply(data_down, function(x) { x[row.names(x) %in% common_names,] })
    data_down <- lapply(data_down, function(x) x[!is.na(x$log2FoldChange),])
    
    #we order stage specifc markers regarding LFC mean of the gene in the different stages in a decrasing way
    data <- seq_len(nrow(data.frame(data_up[1])))
    names(data) <- rownames(data.frame(data_up[1]))
    
    for (g in names(data)){
      data[g] <- mean(unlist(lapply(data_up,function(x) abs(x[g,"log2FoldChange"]))))
      
    }
    data <- sort(data, decreasing = T)
    data <- names(data)
    
    
    
  }else if(updown_markers == "all"){
    #overexpressed genes and infraexpressed are obtainded
    data_up <-getGroupSpecificGenes(select_genes, query_group, updown_markers = "up")
    data_down <- getGroupSpecificGenes(select_genes, query_group, updown_markers = "down")
    data <- as.character(c(data_up,data_down))
    
  }
  
  return(data)
}

# getFilteredDDS <- function(tximport, samples_group, samples_df, filter_min_reads = 5, filter_min_samples  = 1, verbose = TRUE){
#   now.msg(" Filtering DESeqDataSet...", verbose = verbose)
#   
#   # Generation of a factor to contrast
#   cont<-samples_df[,samples_group]
#   df<- cbind(samples_df, design = as.character(cont), stringsAsFactors=FALSE)
#   
#   #DESFromTXI
#   dds <- DESeqDataSetFromTximport(txi = tximport, colData = df, 
#                                   design = ~ design)
#   #estimateSizeFactor
#   dds <- estimateSizeFactors(dds)
#   
#   # Filtering data
#   msg("Fitering data...", verbose = verbose)
#   keep <- rowSums(counts(dds) >= filter_min_reads) >= filter_min_samples
#   filtered_dds <- dds[keep,]
#   
#   now.msg(" DESeqDataSet Filtered", verbose = verbose)
#   return(filtered_dds)
#   
# }

runDESeq2 <- function(sample_names, ref_group, filtered_dds, samples_df, sname_variable, type_lfcShrink = "apeglm", maxpvalue = 0.05, 
                      maxlfc = 0, order_by = "padj", test = "Wald",  verbose = TRUE, ...) {
  
  now.msg(" Starting DESeq...", verbose = verbose)
  
  # Generation of a factor to contrast
  
  cont<-filtered_dds$design
  df <- cbind(samples_df, internal_contrast = as.character(cont), stringsAsFactors=FALSE)
  sample_group_levels <- as.character(df$internal_contrast[samples_df[,sname_variable] %in% sample_names])
  sample_group_levels <- unique(sample_group_levels) 
  
  
  #DESeq
  if(type_lfcShrink == "apeglm"){
    
    msg("Approximate Posterior Estimation for  the  GLM (apeglm) lfcShrink ")
    
    if (test == "LRT"){
      
      msg(" Calculating Likelyhood Ratio Test DESeq ", verbose = verbose)
      
      filtered_dds$design <- droplevels(relevel(filtered_dds$design, ref = ref_group))
      dd <- DESeq(filtered_dds, full = design(filtered_dds), test = "LRT", reduced = ~1, ...)
      
    }else { 
      
      msg(" Calculating Wald Test DESeq ", verbose = verbose)

      filtered_dds$design <- droplevels(relevel(filtered_dds$design, ref = ref_group))
      dd <- DESeq(filtered_dds, full = design(filtered_dds), test = "Wald", ...)
      #dd<- DESeq(filtered_dds, full = design(filtered_dds), test = "Wald")
    }
    
    msg(" Getting DESeq results ", verbose = verbose)
    other_group <- sample_group_levels[which(unique(sample_group_levels != ref_group))]
    res <- results(dd, contrast = c( "design", other_group, ref_group ), alpha = maxpvalue, lfcThreshold = maxlfc)
    res <- lfcShrink(dd, res = res, coef = paste0("design_",other_group, "_vs_", ref_group), type= type_lfcShrink)
    res <- res[order(res[,order_by]),]
    

  }else{
    
    msg("Normaly-distributed prior (Normal) lfcShrink ")
    
    if (test == "LRT"){
      
      msg(" Calculating Likelyhood Ratio Test DESeq ", verbose = verbose)
      
      filtered_dds$design <- droplevels(relevel(filtered_dds$design, ref = ref_group))
      dd <- DESeq(filtered_dds, full = design(filtered_dds), test = "LRT", reduced = ~1, ...)
      
    }else { 
      
      msg(" Calculating Wald Test DESeq ", verbose = verbose)
      
      filtered_dds$design <- droplevels(relevel(filtered_dds$design, ref = ref_group))
      dd <- DESeq(filtered_dds, full = design(filtered_dds), test = "Wald", ...)
      #dd<- DESeq(filtered_dds, full = design(filtered_dds), test = "Wald")
    }
    
    msg(" Getting DESeq results ", verbose = verbose)
    other_group <- sample_group_levels[which(unique(sample_group_levels != ref_group))]
    res <- results(dd, contrast = c( "design", other_group, ref_group ), alpha = maxpvalue, lfcThreshold = maxlfc)
    res <- lfcShrink(dd, res = res, contrast = c( "design", other_group, ref_group ))
    res <- res[order(res[,order_by]),]
    
  }
  
  now.msg("DESeq done", verbose = verbose)
  return (res)
}





# runDESeq2 <- function(sample_names, ref_group, filtered_dds, samples_df, sname_variable, type_lfcShrink = "apeglm", maxpvalue = 0.05, 
#                      maxlfc = 0, order_by = "padj", test = "Wald", verbose = TRUE, ...) {
#   
#   now.msg(" Starting DESeq...", verbose = verbose)
#   
#   # Generation of a factor to contrast
#   cont<-filtered_dds$design
#   df <- cbind(samples_df, internal_contrast = as.character(cont), stringsAsFactors=FALSE)
#   sample_group_levels <- as.character(df$internal_contrast[samples_df[,sname_variable] %in% sample_names])
#   sample_group_levels <- unique(sample_group_levels) 
#    #DESeq
#   
#   if (test == "LRT"){
#     msg(" Calculating Likelyhood Ratio Test DESeq ", verbose = verbose)
#     dd <- DESeq(filtered_dds, full = design(filtered_dds), test = "LRT", reduced = ~1, ...)
#   } else { 
#     msg(" Calculating Wald Test DESeq ", verbose = verbose)
#     dd <- DESeq(filtered_dds, full = design(filtered_dds), test = "Wald", ...)
#     #dd<- DESeq(filtered_dds, full = design(filtered_dds), test = "Wald")
#   }
# 
#   msg(" Getting DESeq results ", verbose = verbose)
#   other_group <- sample_group_levels[which(unique(sample_group_levels != ref_group))]
#   res <- results(dd, contrast = c( "design", other_group, ref_group ), alpha = maxpvalue, lfcThreshold = maxlfc)
#   if(type_lfcShrink == "apeglm"){
#     res <- lfcShrink(dd, res = res, coef = paste0("design_",other_group, "_vs_", ref_group), type= type_lfcShrink)
#   }else{
#     res <- lfcShrink(dd, res = res, c = c( "design", other_group, ref_group ))
#     }
#       
#   res <- res[order(res[,order_by]),]
#   
#   now.msg("DESeq done", verbose = verbose)
#   return (res)
# }







# runDESeq <- function(sample_names, sample_group_levels, ref_group, filtered_dds, samples_df, maxpvalue = 0.05, 
#                      maxlfc = 0, order_by = "padj", filter_min_reads = 1, filter_min_samples  = 1,  
#                      test = "Wald", verbose = TRUE) {
#   now.msg(" Starting DESeq...", verbose = verbose)
#   rnums <- c()
#   for (i in seq_along(sample_names)) {
#     rnums[i] <- which(samples_df$sample.name == sample_names[i])
#   }
#   rnums <- sort(rnums)
#   
#   # Generation of a factor to contrast
#   cont <- rep("unwanted", nrow(samples_df))
#   cont[rnums] <- sample_group_levels
#   samples_df <- cbind(samples_df, internal_contrast = cont, stringsAsFactors=FALSE)
#   
#   #DESFromTXI
#   dds <- DESeqDataSetFromTximport(txi = tximport, colData = samples_df, 
#                                   design = ~ internal_contrast)
#   
#   #estimationSizeFactor
#   dds <- estimateSizeFactors(dds)
#   
#   # data subset
#   dds <- dds[,dds$internal_contrast %in% sample_groups]
#   dds$internal_contrast <- droplevels(relevel(dds$internal_contrast, ref = ref_group))
#   
#   # Filtering data
#   msg("Fitering data...", verbose = verbose)
#   keep <- rowSums(counts(dds) >= filter_min_reads) >= filter_min_samples
#   dds <- dds[keep,]
#   
#   #DESeq
#   
#   if (test == "LRT"){
#     msg(" Calculating Likelyhood Ratio Test DESeq ", verbose = verbose)
#     dds <- DESeq(dds, full = design(dds), test = "LRT", reduced = ~1)
#     msg(" Getting DESeq results ", verbose = verbose)
#     other_group <- sample_group_levels[which(unique(sample_group_levels != ref_group))]
#     res <- results(dds, contrast = c( "internal_contrast", other_group, ref_group ), alpha = maxpvalue, lfcThreshold = maxlfc)
#     res <- lfcShrink(dds, res = res, contrast = c( "internal_contrast", other_group, ref_group ))
#   } else { 
#     msg(" Calculating Wald Test DESeq ", verbose = verbose)
#     dds <- DESeq(dds, full = design(dds), test = "Wald")
#     msg(" Getting DESeq results ", verbose = verbose)
#     other_group <- sample_groups[which(unique(sample_groups != ref_group))]
#     res <- results(dds, contrast = c( "internal_contrast", other_group, ref_group ), alpha = maxpvalue, lfcThreshold = maxlfc)
#     res <- lfcShrink(dds, res = res, contrast = c( "internal_contrast", other_group, ref_group ))
#   }
#       
#   res <- res[order(res[,order_by]),]
#   
#   now.msg("DESeq done", verbose = verbose)
#   return (res)
# }





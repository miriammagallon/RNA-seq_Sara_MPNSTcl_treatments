##################
# Plot Functions #
##################

#########################
# Quality Control plots #
#########################


# MA-plot:

# timepoint_df_position: vector which contains the position of the timepoints to plot. 
# ref_group: the reference timepoint of which the plot is made.
# pseudocounts: rlog(dds) counts it can be filtered pseudocounts or just pseudocounts from rawcounts
# but it is preferably to use filtered pseudocounts in order to compare the data from the same group.


#' MAplot
#' 
#' @description
#' 
#' MAplot function provides a useful overview for an experiment with a two-group comparison
#' 
#' @details
#' 
#' MAplot function plots in the x axis the log-average (A-values) against the Log-Fold Change, as M-values. 
#' M- values  are the count ratio level for each gene between two samples.  
#' A-values are the average count for each gene between two samples.
#' This kind of graphics are useful to visualize if the replicates of samples behave in the same way.
#' Moreover, we see if the data has to be normalized. 
#' In order to do this plot, the data must be filtered.
#' 
#' @usage
#' 
#' MAplot ( pseudocounts, data_df_position, base_filename, base_directory, save = TRUE)
#' 
#' @param pseudocounts (data frame) transformed count data (pseudocounts)
#' @param timepoint_df_position (numeric vector) the factor level position in the data frame.
#' @param base_filename (character) plots name to save
#' @param base_directory(character) the directory where to save all plots generated
#' @param save (logical) whether to save the plots obtained or not(default save = TRUE)
#' 
#' @return
#' 
#' MAplot's main output are png files saved at base_directory or just printed if save = FALSE
#' 
#' @export MAplot


MAplot <- function( pseudocounts, data_df_position, base_filename){
    sampname_time <- colnames(pseudocounts[,data_df_position])
    for( i in 1:(length(data_df_position)-1)){
      nsample <- data_df_position[i]
      for( j in (i+1):length(data_df_position)){
        nsample2 <- data_df_position[j]
        p <- nsample2 +1
        message(nsample, ", ", nsample2)
        M <- pseudocounts[,nsample] - pseudocounts[,nsample2]
        A <- (pseudocounts[,nsample] + pseudocounts[,nsample2])/2 
        df <- data.frame(A, M)
        message (head(df))
        myplot<-print(ggplot(df, aes(x = A, y = M)) + geom_point(size = 1.5, alpha = 1/5) +
                        geom_hline(yintercept =0,color = "blue3") + 
                        stat_smooth(se = FALSE, method = "loess", color = "red3")+ 
                        ggtitle(paste("Time course", base_filename, "days:"), subtitle = paste(sampname_time[i], "vs", sampname_time[j])))
        
        return(myplot)
    }
  } 
}



# Cluster Dendogram

#' getDendogramPlot
#' 
#' @description
#' 
#' getDendogramPlot function provides a herarquical clustering where the relationship between samples is observed.
#' 
#' @details
#' 
#' getDendogramPlot function plots  a herarquical clustering which allows the user to make  an idea
#' about the relation between the studied samples.
#'  
#' @usage
#' 
#' getDendogramPlot(dds, distance_method = "euclidian", plot_title, base_filename, base_directory, save = TRUE, ...)
#' 
#' @param dds (DESeqDataSet) data should be filtered and normalized.
#' @param distance_method (character) method by which data is herarquical clustered. 
#' It could be "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski". (default:"euclidean")
#' @param plot_title (character) the title of the plot. 
#' @param ... further arguments passed to getdendogramPlot
#' 
#' @return
#' 
#' getDendogramPlot's main output is a dendogram plot 
#' 
#' @export getDendogramPlot
#' 

getDendogramPlot <- function(dds, distance_method = "euclidian", plot_title,...){
  
    counts <- counts(dds, normalize = T)
    transposed <- t(counts)
    distance <- dist(transposed, method = distance_method,... )
    clusters <- hclust(distance)
    plot(clusters, main = plot_title, ...)

}

#Heat-mapQC

#' plotHeatmapSamplesQC
#' 
#' @description
#' 
#' plotHeatmapSamplesQC function provides A heatmap of this distance matrix. 
#' 
#' @details
#' 
#' plotHeatmapSamplesQC function provides A heatmap of this distance matrix which
#' gives to the user an overview over similarities and dissimilarities between samples.
#' We have to provide a hierarchical clustering to the heatmap function based on 
#' the sample distances, or else the heatmap function would calculate a clustering based on 
#' the distances between the rows/columns of the distance matrix.
#'  
#' @usage
#' 
#' plotHeatmapSamplesQC(pseudocounts, margins = c(15,15), ...)
#' 
#' @param pseudocounts (data.frame) transformed count data (filtered and normalized)
#' @param margins (numeric)(default: margins = c(15,15))
#' @param ... further arguments passed to plotHeatmapSamplesQC
#' 
#' @return
#' 
#' plotHeatmapSamplesQC's main output are png files saved at base_directory or just printed
#' 
#' @export plotHeatmapSamplesQC
#' 

plotHeatmapSamplesQC <-function(pseudocounts, margins = c(15,15), ...){
  
    mat_dist <- pseudocounts
    mat_dist <- as.matrix(dist(t(mat_dist)))
    mat_dist <- mat_dist/max(mat_dist)
    cim(mat_dist, margins = margins, ...)

}

##################################
#  Differential expression plots #
##################################


# Deg plot MA

#' plotDEGMA
#'
#' @description
#'
#' plotDEGMA function provides a useful overview for an experiment with a two-group comparison
#'
#' @details
#'
#' plotDEGMA function represents each gene with a dot. The x axis is the average expression over all samples,
#' the y axis the log2 fold change between treatment and control. Genes with an adjusted p value below
#' a threshold (by default 0.05) are shown in red.In addition, the most expressed genes can be highlighted.
#' This plot demonstrates that only genes with a large average normalized counts contain sufficient information.
#' to yield a significant call.
#'
#'
#' @usage
#'
#' plotDEGMA(deseq_results, selected_genes, number_genes,...)
#'
#' @param deseq_results (DESeqResults) a DESeqResults object produced by results
#' @param selected_genes (data.frame) Data frame containing the genes already selected by a 0.05 significant level.
#' This data frame is ordered from lowest to highest p-value.
#' @param number_genes (numeric) number of genes that want to be highlighted from the selected_genes data frame.
#' @param title (character) plot title
#' @param ... further arguments passed to plotDEGMA
#'
#' @return
#'
#' plotDEGMA's main output are MA-plots with DEG represented in red and with the most DEG highlighted with a blue circle.
#'
#' @export plotDEGMA
#'


plotDEGMA <- function(deseq_results, selected_genes, number_genes,title, ...){
  
    DESeq2::plotMA(deseq_results, ylim=c(-15,15), xlab = "mean of normalized counts", main= title,...)
    top_gene <- rownames(head(selected_genes, number_genes))
    with (selected_genes[top_gene, ], {
      points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2,...)
      text(baseMean, log2FoldChange, top_gene, pos=2, col="black",cex =0.8, ...)
    })
}


# Heat-map DEG

## blind from rlog function = TRUE should be used for comparing samples in an manner unbiased by prior information on samples.
## bling from rlog function = FALSE should be used for transforming data for downstream analysis, 
## where the full use of the design information should be made. 
## blind=FALSE will skip re-estimation of the dispersion trend, if this has already been calculated. 
## If many of genes have large differences in counts due to the experimental design, it is important to set blind=FALSE for downstream analysis.
## deg_genes: dataframe with all data of selected genes
###Rowv is FALSE this way the order is maintained but, dendogram cant be both sides.

#' getHeatmapPlotDEG
#' 
#' @description
#' 
#' getHeatmapPlotDEG function provides a plot representing and clustering DEG with all samples   
#' 
#' @details
#' 
#' getHeatmapPlotDEG function get a plot representing and clustering DEG with all samples. This function
#' give the oportunity to plot different number of differential expressed genes. In the x axis samples are represented,
#' and in y axis the number of genes that the user like to plot. Blue color  is associated with down regulated genes,
#' while red color is associated with up regulated genes depeneding on LFC associated. The genes are clustered according to
#' its expression.
#' 
#' @usage
#' 
#' getHeatmapPlotDEG(dds, deg_genes, num_genes, margins = c(20, 20),...)
#' 
#' @param dds (DESeqDataSet) a DESeqDataSet object produced by DESeq function
#' @param deg_genes (data.frame) Data frame containing the genes already selected by a significant level. 
#' This data frame is ordered from lowest to highest p-value. 
#' @param num_genes (numeric) number of genes that want to be highlighted from the selected_genes data frame.
#' @param margins (numeric) plot margins (default margins = c(20,20))
#' @param ... further arguments passed to heatMapPlot
#' 
#' @return
#' 
#' getHeatmapPlotDEG's main output are heatmap plots
#' 
#' @export getHeatmapPlotDEG

getHeatmapPlotDEG <- function (dds, deg_genes, num_genes,margins = c(20, 20),...){
  #ds <- rlog(dds)
  ds<-dds
  if(class(ds)== "DESeqDataSet"){
    counts <- data.matrix(assay(ds))
    
  }else if(is.matrix(ds)){
    counts <- ds
    rownames(counts)<- rownames(dds)
  }

  if(is.character(deg_genes)){
    
    heatmap.2(x=counts[rownames(counts)%in% head(deg_genes,num_genes),], col=bluered(75), 
              scale="row", key=TRUE, symkey=FALSE, density.info="none", trace="none", 
               margins = margins, ...)
  }else if (is.data.frame(deg_genes)){
    res <- data.frame(deg_genes)
   
    heatmap.2(x=counts[rownames(head(res, num_genes)), ], col=bluered(75), 
            scale="row", key=TRUE, symkey=FALSE, density.info="none", trace="none", 
            margins = margins, ...)
  }
}

######################################
#   Plot gene Expression functions   #
######################################

#EX:

# A=data.frame("1"=c(2,3,4,5), "2"=c(1,3,5,4), "4"=c(10,3,4,7), "8"=c(5,3,7,9))
# 
# 
# plot(0,0, type= "n", xlim = c(-1,35), ylim= c(0,5000))

# MATRIX

# given one gene we generate a matrix which contains all the information of that gene from
# dds  object (class:DESeq) B = DESeq object (dds)

# dds <- filtered_dds
# points_col = c("red", "blue", "black", "green")
# contrast

#' generationGeneMatrix
#' 
#' @description
#' 
#' generationGeneMatrix function gives a matrix with gene expression information in different samples
#' 
#' @details
#' 
#' generationGeneMatrix function generates a marix were rownames are the sample names and in the colums
#' we find the different points to contrast. Every column will contain at least one count number. This is
#' the first step to get a representation of specific genes expression.
#'  
#' @usage
#' 
#' generationGeneMatrix(dds, gene,  sname_variable, contrast_group, normalize = TRUE)
#' 
#' @param dds (DESeqDataSet) DESeq object with the data filtered 
#' @param gene (character) name of the gene to get its count data. The gene name must be the same gene id as
#' rownames in DESeq object. 
#' @param sname_variable (character) name of the variable of the colData which contains the sample names.
#' @param contrast_group (character) name of the colData variable which contains the contrast points 
#' of expression. 
#' @param normalize (logical) if normalized = TRUE count data will be normalized (default: normalize = TRUE)
#' 
#' @return
#' 
#' generationGeneMatrix's main output is a matrix containing gene expression information of one spcific gene
#' in diferent samples 
#' 
#' @export generationGeneMatrix
#' 

generationGeneMatrix <-function(dds, gene, sname_variable, contrast_group, normalize = TRUE){
  # Getting gene_data
  data <- counts(dds, normalized = normalize)
  gene_data <- data.frame(data[gene,])
  
  # obtaining timepoints
  time <- colData(dds)[,contrast_group]
  col_names <- sort(unique(time))
  
  # Getting the sample names of data
  samples_names <- as.character(colData(dds)[, sname_variable]) 
  row_names <- unique(samples_names)
  
  #Building the matrix from which data will be taken
  ma <- data.frame(matrix(data = NA, nrow = length(row_names), ncol = length(col_names)))
  names(ma)<- col_names
  rownames(ma) <- row_names
  for(i in 1:nrow(gene_data)){
    dd <-colData(dds)[rownames(gene_data)[i],]
    ma[as.character(dd[,sname_variable]), as.character(dd[,contrast_group])] <- gene_data[i,]
    
  }
  return(ma)
}


# A: The matrix got from generationGeneMatrix. It contains the counts in the different contrast_groups. 

# plotLimsAndLabels: Empty plot containing xlim and ylim, as well as the plot title, and xlab and ylab

#' plotLimsAndLabels
#' 
#' @description
#' 
#' plotLimsAndLabels function gives an empty plot containing xlim and ylim inofmation.with gene expression information in different samples
#' 
#' @details
#' 
#' plotLimsAndLabels function generates an empty plot containing xlim and ylim information
#' as well as the title of the graphic. This plot is generated in order to call the function
#' plotPointsAndMeans which will contain the gene expression information an will plot it.
#'  
#' @usage
#' 
#' plotLimsAndLabels(gene_matrix, gene, xlab = "Days", ylab = "Expression level(TPM)" , title = "Expression", ...)
#' 
#' @param gene_matrix (matrix) a matrix generated previusly by generationGeneMatrix function. 
#' @param gene (character) name of the gene studied.
#' @param xlab (character) default "Days"
#' @param ylab (character) default "Expression level (TPM)"
#' @param title (character) plot title (default "Expression")
#' @param ... further arguments passed to plotLimsAndLabels
#' 
#' @return
#' 
#' plotLimsAndLabels's main output is to generate an empty plot containing xlim and ylim information
#' as well as the title of the graphic.
#' 
#' @export plotLimsAndLabels
#'

plotLimsAndLabels <- function(gene_matrix, gene, xlab = "Days", ylab = "Expression level(TPM)" ,  ...){
  plot(0,0, type= "n", xlim = c(min(as.numeric(colnames(gene_matrix))), max(as.numeric(colnames(gene_matrix)))), ylim = c(0,max((gene_matrix),na.rm = TRUE)), xlab = xlab, ylab = ylab , main = paste(gene, "Expression"), ...)
}


# plotsPointsAndMeans: Once we have our empty plot, we proceed to plot the information that we want to.
# (the mean, each sample count value of each group, and a line conecting all mean lines to visualize
# the gene expression profile.)

# plot(0,0, type= "n", xlim = c(min(as.numeric(colnames(A))), max(as.numeric(colnames(A)))), ylim = c(0,max((A),na.rm = TRUE)), xlab = "Days", ylab = "Expression level(TPM)", main = title)
# points_col = as.factor(colData(filtered_dds)$Exp)

#' plotPointsAndMeans
#'
#' @description
#'
#' plotPointsAndMeans function  fill up with information, the empty plot generated by
#' plotLimsAndLables function.
#'
#' @details
#'
#' plotPointsAndMeans function plots lines and points representing the specific gene expression
#' of a query gene.
#'
#' @usage
#'
#' plotPointsAndMeans(gene_matrix, legend_names, legend_title, legend_color,legend_locus, legend_size = 0.5,
#'   points_col = "black", means_col = "black", line_col = "black",
#'   plot_points= TRUE,  plot_means=TRUE, plot_line=TRUE,  means_length = 0.25,
#'   pch = 16, col = "black", plot_legend = TRUE, ...)
#'
#' @param gene_matrix (matrix) a matrix generated previusly by generationGeneMatrix function.
#' @param legend_names (character vector) name of the sample point colours.
#' @param legend_title (character) title of the legend
#' @param legend_color (character) colours of sample points represented in the plot
#' @param legend_size (numeric) size of the legend.
#' @param legend_locus (character) where the user prefers to put the lengend in the plot
#' ("topleft", "topright","bottomleft","bottomright"...)
#' @param points_col (character) point colour(default: points_col = "black")
#' @param means_col (character) mean lines colour (default: means_col = "black")
#' @param line_col (character) colour of the line that join the different means obtained from
#' the query gene expression. (default: line_col = "black")
#' @param plot_points (logical) whether plotting points or not (default: plot_points = TRUE)
#' @param plot_means (logical) whether plotting mean lines or not (default: plot_means = TRUE)
#' @param plot_line (logical) whether plotting line conecting means or not (default: plot_line = TRUE)
#' @param plot_legend  (logical) whether plotting the legend or not (default: plot_legend = TRUE)
#' @param means_length (numeric) length of mean lines. (default: means_length = 0.25)
#' @param pch (nuemric) default pch = 16
#' @param col (character) color by default col = "black"
#' @param ... further arguments passed to plotPointsAndMeans
#'
#' @return
#'
#' plotPointsAndMeans's main output is the representation of points means and lines related to gene expression
#' in different groups.
#'
#' @export plotPointsAndMeans
#'

plotPointsAndMeans <- function(gene_matrix, legend_names, legend_title, legend_color, legend_size = 0.5, legend_locus = "topright", points_col = "black", means_col = "black", line_col = "black", plot_points= TRUE,  plot_means=TRUE, plot_line=TRUE,  means_length = 0.25, pch = 16, points_size = 0.5, col = "black", plot_legend = TRUE, ...){
  timepoints <- as.numeric(gsub("X","",names(gene_matrix)))
  if(plot_points == TRUE){
    for(i in 1:ncol(gene_matrix)){
      points( x =rep(timepoints[i],length(gene_matrix[,i])) , y=gene_matrix[,i], col = points_col, pch = pch, cex = points_size)
    }
  }
  m <- apply(gene_matrix, 2, mean, na.rm = TRUE)
  if(plot_means == TRUE) segments(x0 = timepoints + means_length, x1 = timepoints - means_length, y0 = m, y1 = m, col = means_col, ...)
  if(plot_line == TRUE) lines(x = as.numeric(gsub("X","",names(gene_matrix))), y = m, col = line_col, ...)
  if(plot_legend == TRUE) legend(legend_locus, legend = legend_names , col= legend_color, pch=16, cex = legend_size, title= legend_title, text.font=4, bg= NULL, ...)
}


#plotBarplotGeneExp
plotBarplotGeneExp <- function(dds, gene, sname_variable, contrast_group, bar_order, bar_col, cex.names = 2, cex.axis = 2, ... ){
  A <- generationGeneMatrix(dds = dds, gene = gene, sname_variable = sname_variable, contrast_group = contrast_group, normalize = F)
  A<- A[,bar_order]
  A<- data.matrix(A)
  m <- apply(A, 2, mean, na.rm = TRUE)
  sds <- apply(A,2,sd,na.rm=T)
  bar <- barplot(m, ylim = c(0,max(( m + 2.5*sds ),na.rm = TRUE)), beside = T, main = gene, col =  bar_col, cex.names = cex.names, cex.axis = cex.axis, ...)
  arrows(bar, m-sds,bar, m+sds, angle=90, code=3)
}





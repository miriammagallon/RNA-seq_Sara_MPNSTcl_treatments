#######################################
#         RNA-seq pipeline            #
#######################################
####### Packages Needed
library(DESeq2)
library(tximport)
library(org.Hs.eg.db)
library(yaml)
library(ggplot2)
library(ggbeeswarm)
library(pheatmap)
library(viridis)
if(!require("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(DOSE)
library(gage)
library(apeglm)
library(biomaRt)
library(ggrepel)
library(RColorBrewer)
library(gplots)
library(plyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(regioneR)
library(ggvenn)
library(extrafont)
font_import()
loadfonts("win")
#######################  Loading Functions #########################
source(file= "./RNA-seq_Sara_MPNSTcl_treatments/rna_seq_Functions.R")

#PCAplot deseq function modified
plotPCA <- function (object, intgroup = "condition", pcx1 = "PC1", pcx2="PC2", ntop = 500, returnData = FALSE)
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  pcax <- pca$x
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup,
                                               drop = FALSE])
  if (length(intgroup) > 1) {
    group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }else {
    group <- colData(object)[[intgroup]]
  }
  
  d <- cbind(pcax, group,
             intgroup.df, name = colnames(object))
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:ncol(pcax)]
    return(d)
  }
  ggplot(data = d, aes_string(x = pcx1, y = pcx2, color = "group")) +
    geom_point(size = 3) + xlab(paste0(pcx1,": ", round(percentVar[1] *
                                                          100), "% variance")) + ylab(paste0(pcx2,": ", round(percentVar[2] *
                                                                                                                100), "% variance")) + coord_fixed()
}

#color ramp
map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

#######################     Parameters     ##########################

#### loading the sample data information
sample.data <- read.table(file = "./Sample.info.csv",
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")

sample.group <- "Treatment"
sample.type <- unique(sample.data[,sample.group])



###### Parameters
file.names <- gsub(" ","_",sample.data$Sample.name)
sample.names <- gsub(" ","_",sample.data$Bioinfo.name)


#Salmon parameters
output.suffix <- "_quant"
output.quants <- "./results/Salmon"
# threads <- 8

# Tximport paramenters
orgdb <- org.Hs.eg.db
org.columns <- "SYMBOL"
org.keytype <- "REFSEQ"

# DESeq2 parameters
filt.min.reads <- 5
filt.min.samples <- 1

#Enrichment parameters
goterms <- c("BP", "CC", "MF","KEGG")

color.plate <-  bluered(75)

#Annotation colors
annot.column <- data.frame(Treatment = sample.data[,sample.group], row.names = sample.names)
color.samples <- unique(sample.data$Color.info)
names(color.samples) <- unique(sample.data[,sample.group])
mycolors <- list(Treatment = color.samples)

######### Directories
res.dir <- file.path("results",sample.group)
if(!file.exists(res.dir)) dir.create(res.dir)

pca.dir <-  file.path(res.dir,  "PCAplot")
if(!file.exists(pca.dir)) dir.create(pca.dir)

distance.dir <- file.path(res.dir,  "DistancePlots")
if(!file.exists(distance.dir)) dir.create(distance.dir)

heatmap.QC.dir <- file.path(res.dir, "heatmap.QC")
if(!file.exists(heatmap.QC.dir)) dir.create(heatmap.QC.dir)

heatmap.dir <- file.path(res.dir,"DEG_heatmap")
if(!file.exists(heatmap.dir)) dir.create(heatmap.dir)

DEG.dir <- file.path(res.dir, "DESeqResults")
if(!file.exists(DEG.dir)) dir.create(DEG.dir)

selected.genes.dir <- file.path(res.dir, "DEGgenes")
if(!file.exists(selected.genes.dir)) dir.create(selected.genes.dir)

markers.dir <- file.path(res.dir, "markers")
if(!file.exists(markers.dir)) dir.create(markers.dir)

gsea.dir <- file.path(res.dir,  "GSEAGO")
if(!file.exists(gsea.dir)) dir.create(gsea.dir)
for(term in goterms){
  if(!file.exists(file.path(gsea.dir, term))) dir.create(file.path(gsea.dir, term))
}


annot.dir <- file.path(res.dir, "Annotation")
if(!file.exists(annot.dir)) dir.create(annot.dir)

gene.expression.dir <- file.path(res.dir,  "GeneExpressionPlot")
if(!file.exists(gene.expression.dir)) dir.create(gene.expression.dir)

cluster.dir <- file.path(res.dir,"GeneCluster")
if(!file.exists(cluster.dir)) dir.create(cluster.dir)


##################              Tximport              #####################
salmonquants.fl <-file.path(output.quants, paste0(file.names,"_quant/","quant.sf"))
names(salmonquants.fl)<- file.names
salmonquants.fl <- salmonquants.fl[file.exists(salmonquants.fl)]
txi.salmon <- importQuantsData(quant.files = salmonquants.fl, orgdb = orgdb,
                               orgdb.keytype = orgdb.keytype, org.columns = org.columns,
                               tximpot.type = "salmon")
#Change Sample name
colnames(txi.salmon$abundance) <- sample.names
colnames(txi.salmon$counts) <- sample.names
colnames(txi.salmon$length) <- sample.names
# txi.salmon <- selectDataFromTximport(tximport = txi.salmon, sample_names = sample.names[-c(12,15)])

##################          Quality Control       ######################
sample.data <-sample.data[file.names %in% names(salmonquants.fl),]
# We import the data from tximport output to analyse it using DESeq2
dds.txi <- DESeqDataSetFromTximport(txi = txi.salmon, colData = sample.data, formula(paste0("~", sample.group)))
dim(dds.txi)

# Estimation of Size Factors
dds <- estimateSizeFactors(dds.txi)

# Filtering rawdata: removing data without expression
keep <- rowSums(counts(dds)>=filt.min.reads)>=filt.min.samples
dds <- dds[keep,]
#Normalization
dds.deg <- DESeq(object = dds, test = "Wald")
# dds.degs <- dds.deg[rownames(dds.deg)%in%tf.genes,]
dds.rlog <- rlog(dds.deg)

# Taking a look after filtering and normalizing
# boxplot(counts(dds.deg), main = "BoxPlot_AllSamples_Filtered")

# heatmap QC
mat.dist <- assay(dds.rlog)

mat.dist <- as.matrix(dist(t(mat.dist), method = "euclidean"))
mat.dist <- mat.dist/max(mat.dist)

png(file.path(heatmap.QC.dir, paste0( "QCHeatMap_Euclidean_dist2_AllGenes_ALLSamples_ColumnAnnot_",sample.group,"_validation.png")), width =1000, height = 1000)
pheatmap(mat.dist, annotation_col = annot.column, annotation_row = annot.column, annotation_colors = mycolors, color = viridis(10), fontsize = 10, legend = TRUE, treeheight_row = 50, treeheight_col = 50, lwd = 2, border_color = NA)
dev.off()

#### PCA
# pca.pseudocounts <-rlog(dds.deg)
pca.pseudocounts <- dds.rlog
# colData(pca.pseudocounts)[,sample.group] <- factor(colData(pca.pseudocounts)[,sample.group], levels = names(color.samples))
colData(pca.pseudocounts)[,"group"] <- factor(colData(pca.pseudocounts)[,"PRC2"], levels = c("Mut", "WT"))
colData(pca.pseudocounts)[,"Treatment"] <- factor(x=c("Control","I","A+I","Control","I","A+I","Control","I","Control","I","Control","I","Control","I"), levels= c("Control","I","A+I"))
names(color.samples) <- c("Control","I","A+I")

pcaData <- plotPCA(pca.pseudocounts, intgroup=c("Treatment"), returnData=TRUE)
pcaData$group <- pca.pseudocounts$group
percentVar <- round(100 * attr(pcaData, "percentVar"))

png(filename=file.path(pca.dir, paste0("PCA_PC1-PC2_",sample.group,"_nonames.png")), width = 2000, height = 1500)
svg(filename=file.path(pca.dir, paste0("PCA_PC1-PC2_",sample.group,"_nonames.svg")), width = 40, height = 35)

ggplot(pcaData, aes(PC1, PC2, color = Treatment ))+
  geom_point(size = 44,aes( stroke =3, shape = group, color = Treatment))+ 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  scale_color_manual(name = "Treatment", values = color.samples)+
  scale_shape_manual(name = "PRC2 status", values = c(18,16))+
  # ggrepel::geom_label_repel(aes(label = as.character(rownames(pcaData))),size =10, show.legend = FALSE)+
  theme(legend.title = element_text(size=60), 
        # legend.position = "topright",
        legend.text = element_text( size=40, family = "Arial"), 
        legend.key = element_rect(fill = "white"),
        axis.text=element_text(size=70,family = "Arial"),
        axis.title=element_text(size=70, family = "Arial"),
        plot.title = element_text(size=50, family = "Arial"),
        panel.background = element_rect(fill="white"), 
        panel.grid.major = element_line(colour = "grey50"), plot.margin =margin(30,100,30,30))+
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 30)), 
         shape = guide_legend( override.aes = list(size =30, color ="#7C8e8F", fill = "#7C8e8F" )))


dev.off()  

#################     Differential gene expression      ########################
#Filtering the data
sample.sel.data <- sample.data[sample.data$PRC2 == "Mut",]
# sample.sel.data$Treatment <- factor(sample.sel.data$Treatment,levels=c("Control"))
# sample.sel.data$PRC2 <- factor(sample.sel.data$PRC2, levels=c("Mut"))
# sample.group <- "PRC2"
sample.type <- unique(as.character(sample.sel.data[,sample.group]))
txi.sel <- selectDataFromTximport(tximport = txi.salmon,sample_names = sample.sel.data$Bioinfo.name )
filtered.dds <- getFilteredDDS(tximport = txi.sel, samples_group = sample.group, samples_df = sample.sel.data, filter_min_reads = filt.min.reads, filter_min_samples = filt.min.samples)
dds.rlogsel <-rlog(filtered.dds)

# sample.type <- as.character(unique(sample.sel.data[,sample.group]))
msg(paste("Starting DESeq of", sample.group, "contrast"))
# sample.type <- unique(sample.sel.data$Sample.type)
#Getting the Deseq results
deseq.results <- list()
for(i in seq_len(length(sample.type)-1)) {
  # for(i in seq_len(length(con1))) {
  ref <- sample.type[i]
  # ref <- con1[i]
  
  deseq.results[[ref]] <- list()
  
  for(j in (i+1):length(sample.type)) {
    # for(j in 1:length(con2)) {
    other_group <- sample.type[j]
    # other_group <- con2[j]
    
    # FIPS samples and time %in% c(ref, otro)
    t <- c(ref, other_group)
    deseq_samp <- sample.sel.data[sample.sel.data[,sample.group] %in% t, c("Sample.name", sample.group)]
    
    #runDESeq function
    
    message("RunDeSEQ ", other_group, " vs ", ref)
    
    deseq.results[[ref]][[other_group]] <-runDESeq( sample_names = deseq_samp$Sample.name, sname_variable = "Sample.name", ref_group = ref, samples_df = sample.sel.data, filtered_dds = filtered.dds, type_lfcShrink = "apeglm", test = "Wald")
    
    msg("runDESeq done")
    
  }
  
}

##### DEG ######
selected_genes <- list()

for( i in 1:length(deseq.results)){
  lists <- deseq.results[[i]]
  ref <- names(deseq.results[i])
  contrast <- data.frame()
  for (j in 1:length(lists)){
    contrast <- lists[[j]]
    other_group <- names(lists[j])
    
    # We select those genes which have a pvale <= 0.05 from all analised genes.
    selected_genes[[ref]][[other_group]]<- getTopGenesLFC(deseq_results = contrast, padj_value = 0.05)
    
    selected_genes[[ref]][[other_group]]$genes <- rownames(selected_genes[[ref]][[other_group]])
    # We save the results obtained
    write.table(selected_genes[[ref]][[other_group]],
                file.path(selected.genes.dir, paste0("ResultTable_", ref,"vs", other_group,"_",sample.group,"_Control.csv")),
                sep = "\t", dec = ".", quote = FALSE, col.names = TRUE, row.names = FALSE)
    
    
    write.table(selected_genes[[ref]][[other_group]]$genes,
                file.path(selected.genes.dir, paste0("DEGenes_", ref,"vs", other_group,"_",sample.group,"_Control.csv")),
                sep = "\t", col.names = FALSE, row.names = FALSE)
    
    #selecting those genes with a lfc >2
   genes.lfc2<- rbind(selected_genes[[ref]][[other_group]][selected_genes[[ref]][[other_group]]$log2FoldChange>=2,],selected_genes[[ref]][[other_group]][selected_genes[[ref]][[other_group]]$log2FoldChange<=(-2),])
    
    # We save the results obtained
    write.table(genes.lfc2,
                file.path(selected.genes.dir, paste0("ResultTable_", ref,"vs", other_group,"_",sample.group,"_Control_lfc2.csv")),
                sep = "\t", dec = ".", quote = FALSE, col.names = TRUE, row.names = FALSE)
    
    
    write.table(genes.lfc2$genes,
                file.path(selected.genes.dir, paste0("DEGenes_", ref,"vs", other_group,"_",sample.group,"_Control_lfc2.csv")),
                sep = "\t", col.names = FALSE, row.names = FALSE)
    
    
  }
}


############ Enrichment ##############

genes <- selected_genes$Control$IBET$genes
enrich.terms <- getAllEnrichTermList(genes)

for(i in seq_len(length(goterms))){
  term <- goterms[i]
  if(nrow(data.frame(enrich.terms[[term]]))==0)next
  
  write.table(data.frame(enrich.terms[[term]]), file.path(gsea.dir, term, paste0("Enrichment_", term, "_DEG_ControlvsCombo_",sample.group,".csv")),
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  # write.table(data.frame(enrich.terms[[term]]), file.path(gsea.dir, term, paste0("Enrichment_", term, "_DEG_2KO_MSCvs3KO_MSC_3KO_MSC_lfc2_", sample.group,".csv")),
  # col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  # png(filename = file.path(gsea.dir,term,paste0("Enrich_",term,"DEG_1KO_NC_lfc2", "_1KOvs2KO_NC_",".png")), width = 1000, height = 1000)
  barplot(enrich.terms[[term]],showCategory = 15, font.size = 5)
  ggsave(filename = file.path(gsea.dir,term,paste0("Enrich_",term,"_DEG_ControlvsCombo",".png")), width = 1000, height = 1000,units =  "px")
  # dev.off()
  # cnetplot(enrich.terms[[term]],showCategory = 10, color.params =list(foldChange=genes.sort), categorySize = "pvalue")+
  # theme(plot.title = element_text(size = 5))
  # ggsave(filename = file.path(gsea.dir,term,paste0("Enrich_",term,"_DEG_3KOvsMPNSTcl_MPNSTcl_lfc2_","cnetplot.png")), width = 3500, height = 2000,units =  "px")
  
  # heatplot(enrich.terms$BP,foldChange=genes.sort, showCategory = 30)+
  
  # theme(axis.text= element_text(size = 5))
  
  # ggsave(filename = file.path(gsea.dir,term,paste0("Enrich_",term,"_DEG_1KO_NC_lfc2", "_1KOvsWT_NC_","heatplot.png")), width = 3500, height = 2000,units =  "px")
  
}

################ Heatmap representation #################
data_subset_norm <- t(apply(assay(dds.rlog), 1, cal_z_score))
data_subset_norm <- data_subset_norm[rownames(data_subset_norm) %in%genes,]


genes <-genes[genes %in% rownames(data_subset_norm)]
data_subset_norm <- data_subset_norm[genes,]  
colnames(data_subset_norm)
colnames(data_subset_norm) <- c("SP-10_Control","SP-10_I","SP-10_A+I",
                                "NF1-18B_Control", "NF1-18B_I","NF1_18B_A+I",
                                "S462_Control","S462_I","STS-26T_Control","STS-26T_I",
                                "NF1-09_Control","NF1-09_I", "HS-Sch-2_Control","HS-Sch-2_I")
paletteLength <-length(bluered(80))
myBreaks <- c(seq(min(data_subset_norm), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(data_subset_norm)/paletteLength, max(data_subset_norm), length.out=floor(paletteLength/2)))
annot.column$PRC2_status <- sample.data$PRC2
rownames(annot.column) <- colnames(data_subset_norm)
annot.column$Treatment <- c("Control","I","A+I","Control","I","A+I","Control","I","Control","I","Control","I","Control","I")
mycolors$PRC2_status <- c(Mut="#FF9C78",WT="grey")
names(mycolors$Treatment) <- c("Control","I","A+I")

png(filename = file.path(heatmap.dir,paste0("DEG_IBETvsControl_PRC2_Mut_",sample.group,"_AllSamples_unsupervised_col.png")), width = 1000, height = 800)
pheatmap(data_subset_norm,
         # annotation_row = annotation.row,
         annotation_col = annot.column,
         annotation_colors = mycolors,
         cluster_cols = T,
         cluster_rows = T,
         breaks = myBreaks,
         fontsize = 18, fontsize_row = 7,fontfamily = "arial",
         color = bluered(80),
         show_rownames =T,
         annotation_legend = TRUE, 
         legend = TRUE,
         annotation_names_row = F,
         annotation_names_col = F,
         treeheight_row = 50, 
         treeheight_col = 50,
         lwd = 1)
dev.off() 

################ Barplot and Boxplot representation of specific genes #####################

# genes <- c("MAFA","MGAT4A","LRRC10B","CBLN1","TPPP","VWA5B2","ZDHHC11B","RIMS4","SH3GL2")
# genes <- c("NGFR","TPPP","ZDHHC11B", "PRDM1")
genes <- c("TPPP","NGFR", "PRDM1")

sample.sel.data <- sample.data
sample.sel.data$Treatment <- factor(x=c("Control","I","A+I","Control","I","A+I","Control","I","Control","I","Control","I","Control","I"), levels= c("Control","I","A+I"))
# sample.sel.data <- sample.sel.data[!grepl("HS-Sch-2", sample.sel.data$Bioinfo.name),]
stage <- paste0(sample.sel.data$Treatment,"_",sample.sel.data$PRC2)
stage <- factor(stage,levels= c("Control_Mut","I_Mut","A+I_Mut","Control_WT","I_WT"))

#Normalization
dds.deg <- dds.deg[,!grepl("HS-Sch-2", colnames(dds.deg))]
filt.counts <- counts(dds.deg, normalize = TRUE)

#selecting genes 
filt.counts <- filt.counts[rownames(filt.counts) %in% c(genes),]
genes <- rownames(filt.counts)

# Preparing the data for boxplot representation
fct <- t(filt.counts)
fct <- data.frame(fct)
# length(stage)
# nrow(fct)
fct<- cbind(fct, stage = stage[-c(13:14)])
# fct <- fct[fct$stage%in%c("Control_Mut","IBET_Mut","Combo_Mut"),]
# annot.column
color.data <- sample.sel.data$Color.info
names(color.data) <- stage
color.data[names(color.data)=="Control_WT"] <-"grey"
color.data[names(color.data)=="I_WT"] <- karyoploteR::transparent("#87CDF2")
legend.add <- unique(color.data)
names(legend.add) <- unique(names(color.data))


#Computing mean and deviation of the group


m <- matrix(ncol = ncol(fct)-1, nrow=length(unique(fct$stage)))
desv <- matrix(ncol = ncol(fct)-1, nrow=length(unique(fct$stage)))
g =2
for (g in seq_len(length(unique(stage)))){
  gr <- unique(stage)[g]
  f <-fct[fct$stage %in% gr,]
  f <- f[,-ncol(f)]
  # m <- matrix(ncol = length(genes), nrow=length(unique(fct.list[[gr]]$stage)))
  # desv <- matrix(ncol = length(genes), nrow=length(unique(fct.list[[gr]]$stage)))
  #
  
  # rownames(m)[g] <- gr
  colnames(m) <- colnames(f)
  m[g,] <-colMeans(f)
  # rownames(desv)[g] <-gr
  colnames(desv) <- colnames(f)
  desv[g,] <- apply(f,2,sd)
}
rownames(m) <- unique(fct$stage)
rownames(desv) <- unique(fct$stage)



############################   Barplots ##############
matrix.data <- t(m)
desv.dat <- t(desv)

i=1
for(i in seq_len(length(genes))){
  g <- genes[i]
  if(!(g %in% rownames(matrix.data))) next
  png(filename = file.path(file.path(markers.dir, paste0(g,"_DEG_IBETVSCONTROL_Identity.png"))),
      width = 1500, height = 1000)
  # png(filename = file.path(file.path(markers.dir, paste0(g,"_PCRMarkers_Barplotsclineage.png"))),
  #     width = 1000, height = 800)
  
  par(mar = c(10,10,10,10), mgp=c(7,1,0))
  barCenters <-barplot(matrix.data[g,], beside=TRUE,  
                       cex.names=4, 
                       las=2, 
                       axisnames = FALSE,
                       # ylim=c(0,1.5),
                       ylim=c(0,max(matrix.data[g,]+desv.dat[g,])),#UNNormalized expression
                       cex.axis = 3, 
                       cex.lab = 3,
                       main = g,cex.main =4,
                       # ylab = "Normalized Expression",
                       ylab = "DESeq Expression", #UNnormalized expresion to max value of a gene
                       col = color.data[colnames(matrix.data)])
  arrows( barCenters, matrix.data[g,]-desv.dat[g,],
          barCenters, matrix.data[g,]+desv.dat[g,], angle=90, code= 3, lwd =3)
  # mtext(side=2, line=8, "log2(Expression count)", font=2, cex=3)
  legend(x = 0.3,-0.2,box.lwd = "", xpd = TRUE, fill = color.data[colnames(matrix.data)], legend = names(color.data[colnames(matrix.data)]), horiz = T, cex = 1,col = 2)
  dev.off()
  
  
}

i=1
#BOXPLOT
for(i in seq_len(ncol(fct)-1)){
  # png(filename = file.path("/imppc/labs/eslab/mmagallon/Projects/cNTAP/cNTAP_DEG/Results/LuminexGenes/GenesLuminexCounts",
  # paste0("Boxplot_", colnames(fct)[i], "_normalizedCounts_ALL_luminexGenes_GroupSamples.png")),width = 1200, height = 1000)
  # 
  m <-max(max(fct[,i]))
  
  if(m<=999){
    plus <- 50
    
  } else{
    plus <- 500
    
  }
  pos <- round_any(seq(0,m+plus,by = plus), plus, f = ceiling)
  pos.at  <- c(min(pos),mean(pos),max(pos))
  
  png(filename = file.path(markers.dir,
                           paste0("Boxplot_", colnames(fct)[i], "_DEG_IBETVSCONTROL_Identity_NO_HS-sch2.png")),width = 1200, height = 1000)
  
  par( mar =c(7,15,7,20)+0.1, mgp =c(10,1,0), yaxs = "r")
  bp <- boxplot(fct[,i]~fct[,"stage"], 
                lwd=3,
                xlab = "",
                ylab="",
                main = expression(italic(colnames(fct)[i])),
                boxfill=(legend.add),
                frame=FALSE,
                cex.axis = 5,
                cex.main =3,
                cex.lab = 3,
                xaxt="n",yaxt="n", ylim= c(0, max(m)+plus ))
  
  axis(1, at = c(0:15),labels = rep("",16),lwd.ticks = 0, lwd = 3)
  
  
  # if(m>1000) d <--3
  # if(m<1000) d <- -2
  
  
  # pos <- round(seq(0, m+max(fct[,i]) , by = round(m/3,d=0)),d =0)
  # pos <- round(c(min(pos), mean(pos), max(pos)))
  # pos <- round_any(c(0, m/2, m), 500, f = ceiling)
  
  
  axis(2, at = pos.at, las=2, lwd = 3, cex.axis=5, ylim =c(0,max(m)))
  mtext(side=2, text="DESeq norm counts", line=12, cex = 3)
  # axis(2, at = seq(0,m+100,50000),las=2, lwd = 3, cex.axis=5)
  # axis(2, at = pos, las = 2, lwd = 3, cex.axis =5,ylim =c(0,max(pos)) )
  legend("topright", inset=c(-0.4,0),box.lwd = "", xpd = TRUE, fill = color.data[levels(fct$stage)], legend = names(color.data[levels(fct$stage)]), horiz = F, cex = 3,col = 2)
  # legend(x = 8,y=1000,box.lwd = "", xpd = TRUE, fill = color.data[colnames(matrix.data)], legend = names(color.data[colnames(matrix.data)]), horiz = F, cex = 1,col = 2)
  
  dev.off()
}

i=1
dm <- t(data.matrix(fct)[,1:ncol(fct)-1])
for(i in seq_len(nrow(dm))){
  png(filename = file.path(markers.dir,
                           paste0("Barplot_", rownames(dm)[i], "_NDEG_IBETVSCONTROL_Identity_NOWT.png")),width = 1200, height = 1000)
  
  par( mar =c(20,15,7,20)+0.1, mgp =c(10,1,0), yaxs = "r")
  
  barplot(dm[i,],las=2, main = rownames(dm)[i], col = color.data, cex.axis = 3, cex.names = 2, cex.main =3 )
  legend("topright", inset=c(-0.4,0),box.lwd = "", xpd = TRUE, fill = color.data[levels(fct$stage)], legend = names(color.data[levels(fct$stage)]), horiz = F, cex = 3,col = 2)
  dev.off()
}

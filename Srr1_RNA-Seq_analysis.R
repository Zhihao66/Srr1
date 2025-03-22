library(edgeR)
library(limma)
library(readr)
library(tibble)
library(dplyr)

library(Hmisc)
library(corrplot)
library(pheatmap)
library(ComplexHeatmap)

list.files()
rm(list=ls())
###########################
###Specific gene filtering - keep genes which have a minimum of 5 raw read counts in each of the sample groups
full<-read.delim("all_counts.tsv",row.names=1)
head(full)
colnames(full)
d<-full[c(1:6)]
head(d,5:5)
colnames(d)
#Create new column, "max" for the max read count in each row
d$max<-apply(d,1,max)
head(d)
d1<-d[which(d$max>=5),] #keep rows in which max read count across panel is at least 5 or more
head(d1)
#drop the max column, as its not needed anymore
raw_counts<-d1[c(1:6)]
head(raw_counts)

# define a function to apply the filter to each sample
filter_sample <- function(sample_counts) {
  # count the number of replicates with read count >= 5
  num_pass <- sum(sample_counts >= 5)
  # return 1 if at least 2 replicates pass, 0 otherwise
  if (num_pass >= 2) {
    return(1)
  } else {
    return(0)
  }
}

# apply the filter to each sample, using a rolling window of 3 columns
num_samples <- 2
num_replicates <- 3
pass_filter <- rep(0, nrow(raw_counts))
for (i in 1:num_samples) {
  start_col <- (i-1)*num_replicates + 1
  end_col <- start_col + num_replicates - 1
  sample_counts <- raw_counts[, start_col:end_col]
  sample_pass <- apply(sample_counts, 1, filter_sample)
  pass_filter <- pass_filter + sample_pass
}

# create a new column in the raw_counts table with the pass_filter results
raw_counts$pass_filter <- ifelse(pass_filter > 0, 1, 0)
head(raw_counts)

# write the filtered table to a new file
#write.table(raw_counts, file = "raw_counts.tsv", sep = "\t", quote = FALSE) ##USE THIS TO CHECK THE FILTER eg in excel
#ONLY AFTER MANUALLY CHECKING
raw<-raw_counts[which(raw_counts$pass_filter=="1"),]
head(raw)
d<-raw[c(1:6)]
head(d)
nrow(d)
nrow(raw_counts)
### CALCULATION STEPS
dge <- DGEList(counts=d)
length(rownames(d)) 
head(dge)
#perform scale normalization, adjusting for differences 
dge <- calcNormFactors(dge)
head(dge)
###caculate the CPM value
CPM_allsamples <- cpm(dge)
head(CPM_allsamples)
#write.table(CPM_allsamples, file = "CPM_all_samples.tsv", sep = "\t", quote = FALSE)


#group as multiple rep
group=c(rep("1",3),rep("2",3))
design=model.matrix(~0+group)
colnames(design)=gsub("","",colnames(design))
head(design)

#log transformation with voom
#voom
v=voom(dge,design,plot=TRUE,normalize.method="none") #? voom
head(v)

#MDS
cols<-c(rep("red",3),rep("blue",3))

plotMDS(v, labels = colnames(v),main="min5readcounts_in_any_sample",cex=0.75,col=cols)


###############
###for calculation of DEGs
fit<-lmFit(v,design) # fits row-wise linear models (this will need the voom output file and the design matrix)
fit

#define comparisons as contrast matrices, which are basically AvsB for DEG calculation eg WOR
contrast.matrix= makeContrasts(group1-group2,
                               levels=design)

fit2= contrasts.fit(fit, contrast.matrix)
fit2= eBayes(fit2) #empirical Bayes statistics for differential expression

fit2
#######################################
#######################################
######    1=srr_vs_WT_down	    ###################
#######################################
#######################################
# see contrast matrix. in the example, coef=1 means first comparison (group12-group14 above)
complete_table=topTable(fit2, coef=1,number=nrow(v)) #make sure to change "coef" is changed at each time
head(complete_table)
DEG=complete_table[(complete_table[,"logFC"]>=1|complete_table[,"logFC"]<=(-1))&complete_table[,"adj.P.Val"]<0.05,]
str(DEG)
#write.table(complete_table, file="srr1_vs_WT_complate.tsv", append="false", sep = "\t", quote=F)

upreg=complete_table[complete_table[,"logFC"]>=1&complete_table[,"adj.P.Val"]<0.05,]
#write.table(upreg, file="./upreg_srr1_vs_WT_FC2_P0.05.tsv", append="false", sep = "\t", quote=F)
nrow(upreg)
downreg=complete_table[complete_table[,"logFC"]<=-1&complete_table[,"adj.P.Val"]<0.05,]
#write.table(downreg, file="./downreg_srr1_vs_WT_FC2_P0.05.tsv", append="false", sep = "\t", quote=F)
nrow(downreg)
head(upreg)
head(downreg)
########################### down.  lable the genes #############################################
anno_function<-read.delim("annotation/anno_1.txt")
anno_function
downreg_srr1<-column_to_rownames(downreg_srr1_vs_WT_FPKM,var = "ProtID")
expr_scaled<- t(scale(t(downreg_srr1_vs_WT_FPKM)))
down_TF <- c(anno_function$Srr1_down_TF)
down_TF_genemark <- which(rownames(expr_scaled) %in% down_TF)
down_TF_genemark

############################ down choose the gene to show in the heatmap ####################################
down_sc <- c(anno_function$Srr1_down_Secreted_protein)
genemark_down_sc<- which(rownames(expr_scaled) %in% down_sc )
genemark_down_sc
down_labs<-c(rep("None",nrow(expr_scaled)))
down_labs[genemark_down_sc] <- "Secreted_protein"
down_labs[down_TF_genemark]<- "TF"

down_labs

ha_TF = HeatmapAnnotation(TF = as.factor(down_labs),which = "row",
                          col = list(TF = c("TF" = "blue","None"="white","Secreted_protein" = "white")
                          ),show_legend = F)
ha_SP = HeatmapAnnotation(SP = as.factor(down_labs),which = "row",
                          col = list(SP = c("Secreted_protein" = "green","None"="white","TF"="white")
                          ),show_legend = F)

Heatmap(expr_scaled,
        cluster_columns  = F,
        cluster_rows = T,
        show_row_names = F,
        row_names_gp = gpar(fontsize = 4),heatmap_legend_param = list(title = "Scale")
)+(ha_SP)+(ha_TF)

######################################################
###  DEseq2 for sample correlation   ################
######################################################
######################################################

library(pheatmap)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggrepel)


head(d)
################### set the group data
materials=rep("mycilum",6)
condation<-c("Srr1","Srr1","Srr1","WT","WT","WT")
meta.data=data.frame(materials,condation)
rownames(meta.data)=colnames(d)
all((rownames(meta.data))==colnames(d))

########### make the dds data set
dds_srr1=DESeqDataSetFromMatrix(countData = d,
                                     colData = meta.data,
                                     design= ~condation)

###################caculate the size factor for normalization
dds_srr1<- estimateSizeFactors(dds_srr1)
sizeFactors(dds_srr1)
normlized.count=counts(dds_srr1,normalized=T)
write.csv(normlized.count,"normlized.cout.tsv",col.names=T)

############# samples corration and PCA analysis
vsd_protoplast=vst(dds_srr1,blind = T)
vsd_mat_protoplast=assay(vsd_protoplast)
vsd_cor_protoplast=cor(vsd_mat_protoplast)
View(vsd_cor_protoplast)
pheatmap(vsd_cor_protoplast,annotation = select(meta.data,condation),display_numbers=T)

mds_deseq <- plotMDS(normlized.count)
toplot_deseq <- data.frame(Dim1 = mds_deseq$x, 
                           Dim2 = mds_deseq$y, 
                           Group = as.factor(group),
                           labels = rownames(mds_deseq$cmdscale.out))

ggplot(toplot_deseq, aes(Dim1, Dim2, colour = Group, label=labels)) + 
  geom_point() +
  geom_text_repel(aes(colour = Group), size=3) +
  xlab("Deseq_nomelized_1") +
  ylab("Deseq_nomelized_2") +
  labs(title = paste0("CopciAB","_Scale_log2CPM_Total_gene_reads_MDS_plot"))

plotPCA(vsd_protoplast,intgroup="condation")

dds_srr1=DESeq(dds_srr1)

plotDispEsts(dds_srr1)

######### Volcano plot ######### 
library(EnhancedVolcano)
head(complete_table)
EnhancedVolcano(complete_table,
                lab = rownames(complete_table),
                x = 'logFC',
                y = 'adj.P.Val', pCutoff = 0.05,labSize = 0)

############# GO enrichment #########


library(readr)
library(org.CABnewInterGO1.eg.db)


down_FC_2<-rownames(downreg)

length(down_FC_2)

srr_up_enrichment <- enrichGO(gene=down_FC_2,
                              OrgDb=org.CABnewInterGO1.eg.db,
                              keyType="GID",
                              ont="all",   #CC/BP/MF
                              qvalueCutoff = 0.5,
                              pvalueCutoff =0.05,
                              pAdjustMethod="none",
                              minGSSize = 2,
                              pool=F
)
enrichment_table<-as.data.frame(srr_up_enrichment)

#srr_up_enrichment_plot<-dotplot(srr_up_enrichment,showCategory = 10,color="pvalue",split="ONTOLOGY")
#srr_up_enrichment_plot
#ggsave("srr_down_enrichment_plot_GO_FC2.pdf", plot = srr_up_enrichment_plot, device = "pdf",width = 8,height = 10)



########## KEGG enrichment ########## 
library(ggkegg)
library(ggfx)
library(igraph)
library(tidygraph)
library(dplyr)
library(stringr)
library(readr)
library(tibble)
library(dplyr)
library(readxl)
library(ggplot2)
library(topGO)
library(KEGGREST)
library(Rgraphviz)
library(pathview)
library(GO.db)
library(clusterProfiler)
#In KEGG database there is a Copci dataset, the genome is from Okayama (cci), so first need to find the RBH using MMseq
# 1. ID change of the protein
# CC1G.fa need some ID change from the orignal protein fasta file, and also the ID need to be 5 numbers, #00 should in the beginning if no more than 5 numbers.
#mmeseq to find RBH
#  mmseqs easy-rbh  ./protein_T0.fasta  ./CC1G.fa copci_new_OK.tsv ./tmp


# https://yulab-smu.top/biomedical-knowledge-mining-book/index.html
# https://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/ # local
genetable<-read_tsv("srr_dwonDEG.OK.tsv",col_names = F)
geneList<-genetable$X2
kk <- enrichKEGG(gene         = geneList,
                 organism     = 'cci',
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "none",
                 qvalueCutoff=1,
                 minGSSize = 2)
#dotplot(kk)
#kk
view(as.data.frame(kk))
view(kk@result)
write_tsv(as.data.frame(kk),"KEGG.tsv",col_names = T)



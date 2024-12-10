setwd("C:/Users/Admin/Desktop/Circadian Analyzing data/Hung_diurnal")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")
BiocManager::install("clusterProfiler")
BiocManager::install(c("readr", "pheatmap","RColorBrewer","vsn","dplyr","gplots","AnnotationDbi","AnnotationHub","ggplot2"))
BiocManager::install('EnhancedVolcano')
library(DESeq2)
library(ggplot2)
library(readr)
library(pheatmap)
library(RColorBrewer)
library(vsn)
library(dplyr)
library(gplots)
library(AnnotationDbi)
library(AnnotationHub)
library(tibble)
library(EnhancedVolcano)
##Make the mouse annotation Mouse AH116340
ah <- AnnotationHub()
EnsDb.mouse <- query(ah, c("EnsDb", "Mus musculus"))
EnsDb.mouse
EnsDb.mouse <- EnsDb.mouse[["AH116340"]]
hub <- AnnotationHub()
mouse_resources <- query(hub, "Mus musculus")
mouse_resources # View available resources
resource_id <- "AH116340"
mouse_resource <- hub[[resource_id]]
data <- mouse_resource
### CREATE A DICTIONARY OF GENEID FROM ENSEMBL DATABASE
keys <- keys(data, keytype="GENEID")
maps <- select(data, keys=keys, columns=c("GENENAME", "GENEBIOTYPE"),
               keytype="GENEID")
write.csv(maps, file = "mouse.csv", row.names = FALSE)
mouse <- read.csv("mouse.csv", header = TRUE)


##Convert the data from txt to csv
data <- read.table("counts_clean.txt", header = TRUE, sep = "\t")
write.csv(data, file = "counts_clean_R.csv", row.names = FALSE)
countData <- read.csv("counts_clean_R.csv", header = TRUE, row.names = 1)
countData <-as.matrix(countData)
# Extracting numbers from each column name
new_column_names <- gsub("(\\d+)Aligned.sortedByCoord.out.bam", "\\1", colnames(countData))
# Assign new column names to countData
colnames(countData) <- new_column_names


metaData <- read.csv("metadata.csv", header = TRUE,row.names=1)
metaData
metaData$Group <- as.factor(metaData$Group)
##replaceoutlier
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~ Group)
dds
##Filter only gene > 2000
##keep <- rowSums(counts(dds)) >= 2000
##dds <- dds[keep,]
##Deseq
dds <-DESeq(dds)
result_female <- results(dds, contrast = c("Group","female_active","female_inactive"), alpha = 0.05)
result_male <- results(dds, contrast = c("Group","male_active","male_inactive"), alpha = 0.05)
result_active <- results(dds, contrast = c("Group", "female_active", "male_active"), alpha = 0.05)
result_inactive <- results(dds, contrast = c("Group", "female_inactive", "male_active"), alpha = 0.05)
resultsNames(dds)
## Result in female group active / inactive 
resdata_female <- merge(as.data.frame(result_female), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata_female)[1] <- "EnsemblID"
resdata_female$GeneName <- with(resdata_female, mouse$GENENAME[match(EnsemblID, mouse$GENEID)])
resdata_female <- resdata_female %>% relocate(GeneName, .after = EnsemblID)
resdata_female$GeneType <- with(resdata_female, mouse$GENEBIOTYPE[match(EnsemblID, mouse$GENEID)])
resdata_female <- resdata_female %>% relocate(GeneType, .after = GeneName)
##calculate the mean of group in the comparison
female_group <- colnames(dds)[which(dds$Sex == "female " )]
mean_female <- rowMeans(resdata_female[,colnames(resdata_female) %in% female_group])
resdata_female <- add_column(resdata_female, mean_female, .after = 4)
##resdata_female <- resdata_female[,1:10]
head(resdata_female)
write.csv(as.data.frame(resdata_female), 
          file="Female_Active_vs_Inactive_2.csv")

## Result in male group active / inactive 
resdata_male <- merge(as.data.frame(result_male), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata_male)[1] <- "EnsemblID"
resdata_male$GeneName <- with(resdata_male, mouse$GENENAME[match(EnsemblID, mouse$GENEID)])
resdata_male <- resdata_male %>% relocate(GeneName, .after = EnsemblID)
resdata_male$GeneType <- with(resdata_male, mouse$GENEBIOTYPE[match(EnsemblID, mouse$GENEID)])
resdata_male <- resdata_male %>% relocate(GeneType, .after = GeneName)
##calculate the mean of group in the comparison
male_group <- colnames(dds)[which(dds$Sex == "male" )]
mean_male <- rowMeans(resdata_male[,colnames(resdata_male) %in% male_group])
resdata_male <- add_column(resdata_male, mean_male, .after = 4)
resdata_male <- resdata_male[,1:10]
head(resdata_male)
write.csv(as.data.frame(resdata_male), 
          file="Male_Active_vs_Inactive2.csv")

## Result in active group Female / Male 
resdata_active <- merge(as.data.frame(result_active), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata_active)[1] <- "EnsemblID"
resdata_active$GeneName <- with(resdata_active, mouse$GENENAME[match(EnsemblID, mouse$GENEID)])
resdata_active <- resdata_active %>% relocate(GeneName, .after = EnsemblID)
resdata_active$GeneType <- with(resdata_active, mouse$GENEBIOTYPE[match(EnsemblID, mouse$GENEID)])
resdata_active <- resdata_active %>% relocate(GeneType, .after = GeneName)
head(resdata_active)
write.csv(as.data.frame(resdata_active), 
          file="Active_Female_vs_Male.csv")

## Result in inactive group Female / Male
resdata_inactive <- merge(as.data.frame(result_inactive), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata_inactive)[1] <- "EnsemblID"
resdata_inactive$GeneName <- with(resdata_inactive, mouse$GENENAME[match(EnsemblID, mouse$GENEID)])
resdata_inactive <- resdata_inactive %>% relocate(GeneName, .after = EnsemblID)
resdata_inactive$GeneType <- with(resdata_inactive, mouse$GENEBIOTYPE[match(EnsemblID, mouse$GENEID)])
resdata_inactive <- resdata_inactive %>% relocate(GeneType, .after = GeneName)
head(resdata_inactive)
write.csv(as.data.frame(resdata_inactive), 
          file="Inactive_Female_vs_Male.csv")




####################
#####Graph##########
####################
plotMA(result_female, main=paste0('Female_active_vs_inactive'), ylim=c(-7,7))
plotMA(result_male, main=paste0('Male_active_vs_inactive'), ylim=c(-7,7))

# this gives log2(n + 1)
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))
vsd <- vst(dds, blind=FALSE)
meanSdPlot(assay(vsd))
rld <- rlogTransformation(dds, blind = TRUE)
PCA <-plotPCA(rld, intgroup = c("Group"))
PCA + geom_label(aes(label = Group))
#Plot counts for a single gene. Below is the plot for the gene with the lowest p-value:
geneplot <- "ENSMUSG00000023013"
plotCounts(dds, gene="ENSMUSG00000023013", intgroup="Group", pch = 19, main = mouse[which(mouse$GENEID == geneplot),2])
#Volcano plot 
EnhancedVolcano(result_female,
                lab = mouse$GENENAME[match(rownames(result_female), mouse$GENEID)],
                x = 'log2FoldChange',
                y = 'padj',
                title = "Female active vs inactive time",
                pCutoff = 5e-2,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0)
              
#heatmap of samples
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds))
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

#heatmap and sample to sample distance
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- NULL
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#look for sample outliers
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)




########Plot###########
geneplot <- "ENSMUSG00000050069"
plotCounts(dds, gene=geneplot, intgroup="Group", pch = 19, main = mouse[which(mouse$GENEID == geneplot),2])
print(d)
#####Female Volcano plot########
res <- resdata1[order(resdata1$padj),]
results = as.data.frame(dplyr::mutate(as.data.frame(res), sig=ifelse(res$padj<0.05, "FDR<0.05", "Not Sig")), row.names=rownames(res))
head(results)
DEgenes_DESeq <- results[which(abs(results$log2FoldChange) > log2(1.5) & results$padj < 0.05),]

p = ggplot2::ggplot(results, ggplot2::aes(log2FoldChange, -log10(pvalue))) +
  ggplot2::geom_point(ggplot2::aes(col = sig)) +
  ggplot2::scale_color_manual(values = c("red", "black")) +
  ggplot2::ggtitle("Female_active_vs_inactive")
print(p)
p + ggrepel::geom_text_repel(data=results[1:10, ], ggplot2::aes(label=results[1:10,2]))
#####Male Volcano plot########
res1 <- resdata2[order(resdata2$padj),]
results_male = as.data.frame(dplyr::mutate(as.data.frame(res1), sig=ifelse(res1$padj<0.05, "FDR<0.05", "Not Sig")), row.names=rownames(res1))
head(results_male)
DEgenes_DESeq <- results_male[which(abs(results_male$log2FoldChange) > log2(1.5) & results_male$padj < 0.05),]

p = ggplot2::ggplot(results_male, ggplot2::aes(log2FoldChange, -log10(pvalue))) +
  ggplot2::geom_point(ggplot2::aes(col = sig)) +
  ggplot2::scale_color_manual(values = c("red", "black")) +
  ggplot2::ggtitle("Male_active_vs_inactive")
print(p)
p + ggrepel::geom_text_repel(data=results_male[1:10, ], ggplot2::aes(label=results_male[1:10,2]))


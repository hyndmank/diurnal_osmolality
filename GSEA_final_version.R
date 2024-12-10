setwd("C:/Users/Admin/Desktop/Circadian Analyzing data/Hung_diurnal/GSEA_09-16")
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
library(dplyr)
library(DOSE)
library(RColorBrewer)
BiocManager::install("GOplot",force = TRUE)
library(ReactomePA)
# SET THE DESIRED ORGANISM HERE
##Search in here https://bioconductor.org/packages/release/BiocViews.html#___OrgDb
organism = "org.Mm.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
##Dowload the msigdb database of BroadInstitute 
install.packages("msigdbr")
library(msigdbr)
m_df <- msigdbr(species = "Mus musculus")%>%
  dplyr::select(gs_name, ensembl_gene)
msigdbr_show_species()
## "H" : Hallmark gene sets 
## "C1" : Positional genes sets
## "C2" : Curated gene sets 
## "C3" : Regulatory target gene sets 
## "C4" : Computational gene sets 
## "C5" : Ontology gene sets 
## "C6" : Oncogenic signature gene sets 
## "C7" : immunologic signature gene sets 
## "C8" : cell type signature gene sets 

Mouse_C5 <- msigdbr(species = "Mus musculus", category = "C5") %>%
  dplyr::select(gs_name, ensembl_gene)

#############################  Female #######################################
resdata_Female <- read.csv("Female_active_vs_Female_inactive.csv", header = TRUE)
##resdata_Female_Entrez <- clusterProfiler::bitr(geneID = resdata_Female$EnsemblID,
##fromType="ENSEMBL", toType="ENTREZID",
##OrgDb = organism)
##resdata_Female_Entrez <- merge(resdata_Female, resdata_Female_Entrez, by.x = "EnsemblID", by.y = "ENSEMBL")
##diff <- setdiff(resdata_Female$EnsemblID, resdata_Female_Entrez$EnsemblID)
full_list_Female <- resdata_Female[resdata_Female$mean_female > 10,]
resdata_Female <- resdata_Female[resdata_Female$mean_female > 10 & resdata_Female$padj < 0.05 & !is.na(resdata_Female$padj),]

full_list_Female_Entrez <- clusterProfiler::bitr(geneID = full_list_Female$EnsemblID,
                                                 fromType="ENSEMBL", toType="ENTREZID",
                                                 OrgDb = organism)
DEG_Female <- resdata_Female[which(resdata_Female$padj < 0.05),c("EnsemblID","log2FoldChange")]
DEG_Female_sorted <- DEG_Female[order(DEG_Female$log2FoldChange, decreasing = TRUE),"log2FoldChange"]
names(DEG_Female_sorted) <- DEG_Female[order(DEG_Female$log2FoldChange, decreasing = TRUE),"EnsemblID"]
DEG_Female$group <- "higher_in_Female_active"
DEG_Female$group[DEG_Female$log2FoldChange < 0] <- "higher_in_Female_inactive"

DEG_Female_Entrez <- clusterProfiler::bitr(geneID = DEG_Female$EnsemblID,
                                           fromType="ENSEMBL", toType="ENTREZID",
                                           OrgDb = organism)
DEG_Female_Entrez <- merge(DEG_Female, DEG_Female_Entrez, by.x = "EnsemblID", by.y = "ENSEMBL")

############################## Male ########################################
resdata_Male <- read.csv("Male_active_vs_Male_inactive.csv", header = TRUE)
full_list_Male <- resdata_Male[resdata_Male$mean_male > 10,]
resdata_Male <- resdata_Male[resdata_Male$mean_male > 10 & resdata_Male$padj < 0.05 & !is.na(resdata_Male$padj),]
full_list_Male_Entrez <- clusterProfiler::bitr(geneID = full_list_Male$EnsemblID,
                                               fromType="ENSEMBL", toType="ENTREZID",
                                               OrgDb = organism)
DEG_Male <- resdata_Male[which(resdata_Male$padj < 0.05),c("EnsemblID","log2FoldChange")]
DEG_Male_sorted <- DEG_Male[order(DEG_Male$log2FoldChange, decreasing = TRUE),"log2FoldChange"]
names(DEG_Male_sorted) <- DEG_Male[order(DEG_Male$log2FoldChange, decreasing = TRUE),"EnsemblID"]
DEG_Male$group <- "higher_in_Male_active"
DEG_Male$group[DEG_Male$log2FoldChange < 0] <- "higher_in_Male_inactive"

DEG_Male_Entrez <- clusterProfiler::bitr(geneID = DEG_Male$EnsemblID,
                                         fromType="ENSEMBL", toType="ENTREZID",
                                         OrgDb = organism)
DEG_Male_Entrez <- merge(DEG_Male, DEG_Male_Entrez, by.x = "EnsemblID", by.y = "ENSEMBL")


##Run the input from GSEA_09-16.R first
up_female <- DEG_Female[DEG_Female$log2FoldChange > 0, "EnsemblID"]
down_female <- DEG_Female[DEG_Female$log2FoldChange < 0, "EnsemblID"]
up_female_entrez <- DEG_Female_Entrez[DEG_Female_Entrez$log2FoldChange > 0, "ENTREZID"]
down_female_entrez <- DEG_Female_Entrez[DEG_Female_Entrez$log2FoldChange <0, "ENTREZID"]
up_male <- DEG_Male[DEG_Male$log2FoldChange > 0, "EnsemblID"]
down_male <- DEG_Male[DEG_Male$log2FoldChange < 0, "EnsemblID"]
up_male_entrez <- DEG_Male_Entrez[DEG_Male_Entrez$log2FoldChange > 0, "ENTREZID"]
down_male_entrez <- DEG_Male_Entrez[DEG_Male_Entrez$log2FoldChange <0, "ENTREZID"]


list_up <- list(Female=up_female,Male=up_male)
list_down <- list(Female=down_female,Male=down_male)
list_up_entrez <- list(Female=up_female_entrez,Male=up_male_entrez)
list_down_entrez <- list(Female=down_female_entrez,Male=down_male_entrez)
background_gene <- unique(c(full_list_Female$EnsemblID,full_list_Male$EnsemblID))
background_gene_entrez <- unique(c(full_list_Female_Entrez$ENTREZID,full_list_Male_Entrez$ENTREZID))
enrichGO_Up <- compareCluster(  geneClusters = list_up,
                                universe=background_gene,
                                ont           = "ALL",
                                OrgDb = organism,
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05,
                                readable      = TRUE,
                                fun =  enrichGO,
                                keyType = "ENSEMBL")
write.csv(enrichGO_Up, "enrichGO_Up_in_female_male.csv")

cnetplot(enrichGO_Up,showCatergory = 20, cex.params = list(category_node = 3, gene_node = 1, category_label = 2, gene_label =1, 
                                                               color_category = "red", color_gene = "blue"), max.overlaps = 10)

dotplot(enrichGO_Up, showCategory = 6,
        ##split = "ONTOLOGY", 
        label_format = 30) +
  ggtitle("Higher enrichment in active time in Female and Male") +
  theme(plot.title = element_text(hjust = 1))



enrichKEGG_Up <- compareCluster(geneClusters = list_up_entrez,
                                    universe=background_gene_entrez,
                                    ##ont           = "ALL",
                                    ##https://www.genome.jp/kegg/tables/br08606.html
                                    organism = "mmu",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.01,
                                    qvalueCutoff  = 0.05,
                                    ##readable      = TRUE,
                                    fun = enrichKEGG)
#the default is ENTREZID, change the keyType cause error)
enrichKEGG_Up <- setReadable(enrichKEGG_Up, OrgDb = organism, keyType="ENTREZID")   
enrichKEGG_Up@compareClusterResult$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", enrichKEGG_Up@compareClusterResult$Description)
write.csv(enrichKEGG_Up, "enrichKEGG_Up.csv")
cnetplot(enrichKEGG_Up,cex.params = list(category_node = 2, gene_node = 1, category_label = 2, gene_label =1, 
                                             color_category = "red", color_gene = "blue"))
dotplot(enrichKEGG_Up, showCategory = 16)

enrichGO_Down <- compareCluster(geneClusters = list_down,
                                universe=background_gene,
                                ont           = "ALL",
                                OrgDb = organism,
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05,
                                readable      = TRUE,
                                fun =  enrichGO,
                                keyType = "ENSEMBL")
write.csv(enrichGO_Down, "enrichGO_Down_in_female_male.csv")

cnetplot(enrichGO_Down, showCategory = 12,  # Increase the number of categories
         foldChange = NULL,
         cex.params = list(category_node = 3, gene_node = 1, 
                           category_label = 2, gene_label = 1, 
                           color_category = "red", color_gene = "blue"), 
         max.overlaps = 20)

dotplot(enrichGO_Down_simplified, showCategory = 11,
        , label_format = 30) +
  ggtitle("Higher enrichment in inactive time in Female and Male") +
  theme(plot.title = element_text(hjust = 1))

###Simplified######
enrichGO_Down_simplified <- simplify(enrichGO_Down, cutoff=0.7, by="p.adjust", select_fun=min, measure = "Wang")
write.csv(enrichGO_Down_simplified, "enrichGO_Down_simplified_in_female_male.csv")
##################

enrichKEGG_Down <- compareCluster(geneClusters = list_down_entrez,
                                universe=background_gene_entrez,
                                ##ont           = "ALL",
                                ##https://www.genome.jp/kegg/tables/br08606.html
                                organism = "mmu",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05,
                                ##readable      = TRUE,
                                fun = enrichKEGG)
#the default is ENTREZID, change the keyType cause error)
enrichKEGG_Down <- setReadable(enrichKEGG_Down, OrgDb = organism, keyType="ENTREZID")
enrichKEGG_Down@compareClusterResult$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", enrichKEGG_Down@compareClusterResult$Description)

write.csv(enrichKEGG_Down, "enrichKEGG_Down.csv")
cnetplot(enrichKEGG_Down,cex.params = list(category_node = 1.5, gene_node = 0.75, category_label = 2, gene_label =1, 
                                         color_category = "red", color_gene = "blue"))
dotplot(enrichKEGG_Down, showCategory = 16)




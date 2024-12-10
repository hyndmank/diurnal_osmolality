resdata_active_vs_inactive_female <- read.csv("Female_active_vs_Female_inactive.csv", header = TRUE)
ENTREZID  <- clusterProfiler::bitr(geneID = resdata_active_vs_inactive_female$EnsemblID,
                                   fromType="ENSEMBL", toType="ENTREZID",
                                   OrgDb = organism)
resdata_active_vs_inactive_female <- resdata_active_vs_inactive_female %>%
  left_join(ENTREZID, by = c("EnsemblID" = "ENSEMBL"))

full_list_Female <- resdata_active_vs_inactive_female[resdata_active_vs_inactive_female$mean_female > 10,]

DEG_active_vs_inactive_female <- full_list_Female[which(full_list_Female$padj < 0.05 &
                                                          !is.na(full_list_Female$padj))
                                                  ,c("EnsemblID","log2FoldChange","ENTREZID")]

dim(DEG_active_vs_inactive_female)

length(unique(DEG_active_vs_inactive_female$EnsemblID))

DEG_active_vs_inactive_female_sorted <- DEG_active_vs_inactive_female[order(DEG_active_vs_inactive_female$log2FoldChange, 
                                                                            decreasing = TRUE),]
DEG_active_vs_inactive_female_sorted_ENTREZID <- DEG_active_vs_inactive_female_sorted[
  which(!is.na(DEG_active_vs_inactive_female_sorted$ENTREZID)),"log2FoldChange"]
names(DEG_active_vs_inactive_female_sorted_ENTREZID) <- DEG_active_vs_inactive_female_sorted[
  which(!is.na(DEG_active_vs_inactive_female_sorted$ENTREZID)),"ENTREZID"]
DEG_active_vs_inactive_female_sorted_ENSEMBL <- DEG_active_vs_inactive_female_sorted[,"log2FoldChange"]
names(DEG_active_vs_inactive_female_sorted_ENSEMBL) <- DEG_active_vs_inactive_female_sorted[,"EnsemblID"]

length(DEG_active_vs_inactive_female_sorted_ENSEMBL)
###############
kk2 <- gseKEGG(geneList     = DEG_active_vs_inactive_female_sorted_ENTREZID,
               organism     = "mmu",
               #nPerm        = 1000,
               minGSSize    = 10,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               keyType       = "ncbi-geneid")
kk2 <- setReadable(kk2, OrgDb = organism, keyType="ENSEMBL")  
dotplot(kk2, showCategory = 20)
cnetplot(kk2,cex.params = list(category_node = 3, gene_node = 1, category_label = 2, gene_label =1, 
                                   color_category = "red", color_gene = "blue"))
write.csv(kk2, "heat1h_control_female_gseKEGG.csv")




library(pathview)
# Define the KEGG pathway ID
pathway_id <- "mmu04962"  # Mouse "Vasopressin-regulated water reabsorption"
kegg_genes <- bitr_kegg(pathway_id, fromType = "kegg", toType = "symbol", organism = "mmu")
print(kegg_genes)

# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=DEG_active_vs_inactive_female_sorted_ENTREZID, pathway.id="mmu04962", species = "mmu")
print(dme)
# Produce a different plot (PDF) (not displayed here)
dme <- pathview(gene.data=DEG_active_vs_inactive_female_sorted_ENTREZID, pathway.id="mmu04962", species = kegg_organism, kegg.native = F)


##########Vasopressin pathway ############
if (!requireNamespace("KEGGREST", quietly = TRUE)) {
  install.packages("KEGGREST")
}
library(KEGGREST)
pathway_id <- "mmu04962"  # Vasopressin-regulated water reabsorption in mice
pathway_data <- keggGet(pathway_id)
genes <- pathway_data[[1]]$GENE
gene_symbols <- genes[seq(1, length(genes), by = 2)]
length(gene_symbols)
##########################################
results_table <- data.frame(Gene = character(), log2 = numeric(), padj = numeric(), genenam = character(), stringsAsFactors = FALSE)

# Loop through each gene in the vector
for (i in gene_symbols) {
  # Extract the log2FoldChange values for active_fe and HS
  log2 <- resdata_active_vs_inactive_female[resdata_active_vs_inactive_female$ENTREZID == i
                                                 & !is.na(resdata_active_vs_inactive_female$ENTREZID), "log2FoldChange"]
  padj <- resdata_active_vs_inactive_female[resdata_active_vs_inactive_female$ENTREZID == i
                                            & !is.na(resdata_active_vs_inactive_female$ENTREZID), "padj"]
  genenam <- resdata_active_vs_inactive_female[resdata_active_vs_inactive_female$ENTREZID == i
                                               & !is.na(resdata_active_vs_inactive_female$ENTREZID), "GeneName"]
  # Handle cases where no match is found (NA)
  log2 <- ifelse(length(log2) == 0, NA, log2)
  padj <- ifelse(length(padj) == 0, NA, padj)
  genenam <- ifelse(length(genenam) == 0, NA, genenam)
  # Add the result as a new row in the data frame
  results_table <- rbind(results_table, data.frame(Gene = i, log2 = log2, padj = padj, genename = genenam))
}
dim(results_table)

ggplot(results_table, aes(x = genename, y = log2, fill = padj < 0.05)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
  labs(
    title = "Bar Plot of Gene Expression",
    x = "Gene Name",
    y = "Log2 Fold Change",
    fill = "Significant (padj < 0.05)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for clarity

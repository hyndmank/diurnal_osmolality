if (!requireNamespace("KEGGREST", quietly = TRUE)) {
  install.packages("KEGGREST")
}
library(KEGGREST)

# Define the KEGG pathway ID
pathway_id <- "mmu04962"  # Vasopressin-regulated water reabsorption in mice

# Retrieve the pathway details
pathway_data <- keggGet(pathway_id)

# Extract the list of genes
genes <- pathway_data[[1]]$GENE

# Extract gene symbols (genes alternate with descriptions in the list)
gene_symbols <- genes[seq(1, length(genes), by = 2)]

# Print the list of gene symbols
print(gene_symbols)

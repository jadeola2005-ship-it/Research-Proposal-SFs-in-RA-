# Load required packages 
library(tximport)         # import quantification files 
library(DESeq2)           # differential expression 
library(clusterProfiler)  # Pathway enrichment 
library(org.Hs.eg.db)     # human gene annotation 
library(STRINGdb)         # protein - protein interactions 
library(ggplot2)          # Visualisation 

# Import quantification files 
files <- c(
"WT_1/quant.sf","WT_2/quant.sf","WT_3/quant.sf",
"R620W_1/quant.sf","R620W_2/quant.sf","R620W_3/quant.sf",
"KO_1/quant.sf","KO_2/quant.sf","KO_3/quant.sf"
)

txi <- tximport(files, type = "salmon", txOut = FALSE)

# Metadata for sample groups 
coldata <- data.frame(
 row.names = colnames(txi$counts),
 condition = factor(c(rep("WT",3), rep("R620W",3), rep("KO",3)))
)

# Differential expression with DESeq2 
dds <- DESeqDataSetFromTximport(txi, coldata, design = ~ condition)
dds <- DESeq(dds)

# Example contrast: R620W vs WT
res_R620W <- results(dds, contrast = c("condition","R620W","WT"))

# Extract signifcant genes (FDR < 0.05)
sig_genes <- rownames(res_R620W)[res_R620W$padj < 0.05]

# Pathway enrichment (Gene Ontology)
go_results <- enrichGO(
  gene          = sig_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH"
)


# Protein to Protein interactionb network (STRINGdb)
string_db <- STRINGdb$new(version = "12", species = 9606, score_threshold = 400)
mapped <- string_db$map(data.frame(gene = sig_genes), "gene", removeUnmapped = TRUE)
ppi_network <- string_db$get_interactions(mapped$STRING_id)

# PCA plot (Visualising separation between groups)
### PCA Analysis Based on refine.bio Example ###

# Load packages
library(DESeq2)
library(ggplot2)

# Variance stabilising transformation (VST)
vsd <- vst(dds, blind = TRUE)

# Compute PCA using DESeq2 helper function
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

# Calculate percent variance explained (as in refine.bio)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# PCA plot
ggplot(pcaData, aes(PC1, PC2, colour = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of WT, R620W, and KO Synovial Fibroblasts")

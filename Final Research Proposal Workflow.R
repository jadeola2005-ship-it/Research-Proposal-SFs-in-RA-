# Step 1 write a Salmon-format quant.sf file
write_quant <- function(df, path) {
  # Ensure that columns are exactly in the correct order
  df <- df[, c("Name","Length","EffectiveLength","TPM","NumReads")]
  
  # Write tab-separated file with EXACT Salmon formatting
  write.table(df,
              file = path,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE)
}

simulate_quant <- function(path, condition){
  
  # Create 500 gene names
  genes <- paste0("Gene", 1:500)
  
  # Baseline random values
  Length          <- sample(500:2000, 500, TRUE)
  EffectiveLength <- sample(400:1800, 500, TRUE)
  TPM             <- runif(500, 0, 20)
  NumReads        <- rpois(500, 80)
  
  # R620W variant → autophagy + mitochondrial + inflammatory upregulation
  if(condition == "R620W"){
    autophagy <- 1:15       # pretend these are ATG5, BECN1, LC3B
    mito      <- 16:30      # pretend these are SOD2, PINK1
    inflam    <- 31:45      # pretend these are IL6, IL8, CCL2
    
    TPM[autophagy] <- TPM[autophagy] + runif(15, 10, 25)
    TPM[mito]      <- TPM[mito]      + runif(15, 15, 30)
    TPM[inflam]    <- TPM[inflam]    + runif(15, 20, 40)
    
    NumReads[autophagy] <- NumReads[autophagy] + rpois(15, 50)
    NumReads[mito]      <- NumReads[mito]      + rpois(15, 60)
    NumReads[inflam]    <- NumReads[inflam]    + rpois(15, 80)
  }
  
  # KO samples show loss of PTPN22 expression
  if(condition == "KO"){
    TPM[100]      <- runif(1, 0, 0.5)  
    NumReads[100] <- rpois(1, 2)       
  }
  
  # Build quant.sf-like table
  quant <- data.frame(
    Name            = genes,
    Length          = Length,
    EffectiveLength = EffectiveLength,
    TPM             = TPM,
    NumReads        = NumReads,
    stringsAsFactors = FALSE
  )
  
  # Save quant.sf files
  write_quant(quant, file.path(path, "quant.sf"))
}


 # Step 3 created folders and generate simulated data 
sample_dirs <- c(
  "WT_1","WT_2","WT_3",
  "R620W_1","R620W_2","R620W_3",
  "KO_1","KO_2","KO_3"
)

# Create folders
for (folder in sample_dirs) dir.create(folder, showWarnings = FALSE)
 
# Create folders
for (folder in sample_dirs) dir.create(folder, showWarnings = FALSE)

# Generate simulated quant.sf files
simulate_quant("WT_1", "WT")
simulate_quant("WT_2", "WT")
simulate_quant("WT_3", "WT")

simulate_quant("R620W_1", "R620W")
simulate_quant("R620W_2", "R620W")
simulate_quant("R620W_3", "R620W")

simulate_quant("KO_1", "KO")
simulate_quant("KO_2", "KO")
simulate_quant("KO_3", "KO")
  
 

# Step 4 — Build file list for R  
files <- file.path(
  c("WT_1","WT_2","WT_3",
    "R620W_1","R620W_2","R620W_3",
    "KO_1","KO_2","KO_3"),
  "quant.sf"
)

# Check all files exist
file.exists(files)


# Step 5 — Load required libraries   

library(tximport)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(STRINGdb)
library(ggplot2)


# Step 6 — Import simulated RNA-seq    
txi <- tximport(files, type = "salmon", txOut = TRUE)

#confirm structure 
names(txi)

#check dimensions 
dim(txi$counts)


# Step 7 — Build DESeq2 object 
coldata <- data.frame(
  row.names = colnames(txi$counts),
  condition = factor(c(rep("WT",3), rep("R620W",3), rep("KO",3)))
)

dds <- DESeqDataSetFromTximport(txi, coldata, design = ~ condition)
dds <- DESeq(dds)

# R620W vs WT results
res_R620W <- results(dds, contrast = c("condition","R620W","WT"))

# Significant genes
sig_genes <- rownames(res_R620W)[which(res_R620W$padj < 0.05)]


#Step 8 — GO enrichment analysis 
go_results <- enrichGO(
  gene = sig_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH"
)


#Step 9 — Stringdb protein interaction network 
string_db <- STRINGdb$new(version = "12", species = 9606, score_threshold = 400)

mapped <- string_db$map(
  data.frame(gene = sig_genes),
  "gene",
  removeUnmapped = TRUE
)

ppi_network <- string_db$get_interactions(mapped$STRING_id)


# Step 10 — PCA plot of all samples 
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, colour = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of WT, R620W, and KO Synovial Fibroblasts")



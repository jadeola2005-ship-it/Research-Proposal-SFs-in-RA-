Project Overview 
This repositroy contains an exaple RNA smequencing analysis workflow designed to support the research proposal in investigating the role of the PTPN22 R620W variant in synovial fibroblasts (SFs) in Rheumatoid arthritis (RA)
Because no real RNA seq data is availabel for this proposal the worflow uses stimulated sequencing data to demonstrate how the analysis would be performed in a real experimental setting.
The pipeline reflects standard RNA Seq bioinformatics practice and is intened to show and understand the methods, reproducibility and feasibility.

Research aim 
demonstrate how RNA Seq can be used to compare gene expression between:
WT synovial fibrobalsts (normal PTPN22)
PTPN22 R620W variant fibroblasts 
PTPN22 Knockout (KO) fibroblasts 

with a focus on pathways related to:
Autophagy 
Mitochondrial stress 
Inflammatory signalling 

The R script performs the following steps:
Stomulates RNA quantification files 
- Genereate realistic Salmon format quant.sf files
- creates 9 Samples (3WT, 3R620W and 3KO)
R620W samples show increased autophagy, mithchondrial stress and inflammatory signals
KO samples mimic loss of PTPN22 function

Imports data using tximport 
- Read stimulated quant.sf files into R
- Uses tx0ut = TRUE because the data are stimulated and do not contain transcript gene mappings
This steps follows the frameowrk described by (Soneson, 2025)

Performs differential expression analysis with DESeq2 
- Builds a DESeq2 dataset
- compares gene expression between WT,R620W and KO fibroblasts
- identifies genes altered by the PTPN22 R620W variant
This follows methods described by (Love,2014)

Functional pathway analysis (Gene Ontology) 
Uses clusterProfiler to identify enriched biological processes 
Focuses on Pathways relevant to RA pathology 
Similar to methods from (Yu, 2012)

PPI (Protein - protein interaction analysis)
Uses STRINGfb to construct predicted interaction networks 
Highlights potential regulatory hubs altered by the R620W variant 
Similar methods to (Szklarczyk, 2015)

Principle component analysis (PCA)
- Visualises global gene expression difference between groups
- confirms separation between WT, R620W and KO fibroblasts
- Demonstrates that PTPN22 genotype influences transcriptomic state
PCA is used for variance analysis and quality control not correlation testing

Folder Structure 
when the script is run these files are automatically created 
<img width="251" height="368" alt="image" src="https://github.com/user-attachments/assets/c177dca3-41ec-4339-9f6b-9d54e97eb5db" />
Each folder contains a quant.sf file 

The R packages used are:
tximport - import Salmon quantifications 
DESeq2 - Differential expression analysis 
ClusterProfiler - Pathway enrichment 
STRINGdb - Pathway enrichment 
ggplot2 - visualisation 

Key notes 
This is a simulated dataset 
Biological conclusions could be drawn but the output is hypothetical 
The script shows how the analysis should be conducted and the results that could be expected 
tx0ut = TRUE is needed as transcript to gene mappings are not available in simulated data 

Simulated RNA Seq purpose 
- DEmonstrate a complete and realsitic bioinformatics pipeline
- Show how hypothesis would be tested using RNA Seq
- provide reproducible and transpoarent methods
- Prepare the analysis frameowrk for future datasets that are real

Summary 
This repository provides a clear and reproducible example RNA seq worflow that supports the research proposal on the role of PTPN22 R620W in SFs











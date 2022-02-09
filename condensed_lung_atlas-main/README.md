# Condensed atlas for lung development
## Purpose
The purpose of this repo is to create a condensed atlas of scRNA-seq data for the developing lung. 

## Workflow
1. Ingest multiple datasets over a common tissue (e.g. lung), all sharing common metadata fields and in counts per million (CPM). 
   1. Metadata
      1. Dataset*
      2. Cell type (original)
      3. Cell type (harmonized)*
      4. Lineage (endothelial, epithelial, mesenchymal, immune)
      5. Timepoint*
      6. Sex 
      7. Technology 
      8. Individual
      
      "*" Denotes what is currently moved to condensed atlas

2. Use *northstar* to call a cell type for every cell in each dataset and only looking for nearest neighbors in the atlas provided
   1. Atlas is the combination of these three datasets 
      1. [Immune](https://pubmed.ncbi.nlm.nih.gov/32484158/)
      2. [Endothelial](https://www.biorxiv.org/content/10.1101/2021.04.27.441649v1.full)
      3. [Mesenchymal](https://www.biorxiv.org/content/10.1101/2021.05.19.444776v1) 
      4. Epithelial (to be added)
         1. Need to decide on source
3. Concatenate .h5ad of all datasets together after cell typing
4. Create .h5 file for all data/stats of interest from merged .h5ad
   1. For each field of metadata
      1. Cell count
      2. Gene expression matrix
      3. Proportion expression matrix
   2. Number of metadata fields in .h5 = n! 
      1. n = length of metadata column list
      2. Every combination of metadata gets a key
         
      


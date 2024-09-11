### About the tool
##### ScInfeR: a graph-based cell type annotation toolkit for single-cell RNA-seq, ATAC-seq, and spatial omics
##### Brief overview of the tool algorithm
1. ScInfeR can annotate cells using user-defined marker sets, scRNA-seq references, or both to annotate cells in scRNA-seq, scATAC-seq, and spatial omics datasets.  
2. ScInfeR implements two rounds of annotation strategy for cell-type assignment. 
3. First, the tool annotates the cell clusters by correlating the cluster-specific markers with the cell-type-specific markers in the cell-cell similarity graph. These cell-type-specific marker genes can be either user-defined, extracted by ScInfeR from scRNA-seq reference data, or a combination of both. For scRNA-seq as a reference, ScInfeR extracts cell-type markers by considering both the global and local specificity of markers. 
4. In the second round, the tool annotates the subtypes and clusters containing multiple cell types in a hierarchical manner. In this step, the tool uses a framework adapted from the message-passing layer in the graph neural network to annotate each cell individually.
### How to install 
#### using devtools
```
devtools::install_github("swainasish/ScInfeR")
```
#### Or you can load the R script as source
```
source("https://raw.githubusercontent.com/swainasish/ScInfeR/master/R/base.R")
```
\name{predict_celltype_scRNA_seurat}
\alias{predict_celltype_scRNA_seurat}
\title{Cell type prediction using a Seurat scRNA-seq object}
\usage{
predict_celltype_scRNA_seurat(s_object,
                              group_annt,
                              ct_marker_df,
                              subtype_present = F,
                              subtype_info = F,
                              assay_name="RNA",
                              slot_name="counts",
                              reduction="umap",
                              own_weightage = 0.5,
                              n_neighbor=10)
}
\description{
This function predicts celltype annotations from a Seurat scRNA-seq object.
}
\arguments{
\item{s_object}{scRNA-seq Seurat object}
\item{group_annt}{Cluster information obtained from findclusters function in Seurat}
\item{ct_marker_df}{dataframe of celltype markers, should have 3 columns having names: "celltype","marker","weight"}
\item{subtype_present}{Whether the subtypes present in the ct_marker_df}
\item{subtype_info}{if the subtype_present=T, subtype information in a list form.
Example: subtype_stat <- list("T-cell"=c("CD8+ T-cell","CD4+ T-cell"),
              "B cell"=c("Naive B cell","Memory B cell"))}
\item{assay_name}{Assay name in the Seurat object}
\item{slot_name}{Slot name, where the log normalised count matrix present in the Seurat object}
\item{reduction}{Name of the reduction method used during the downstream analysis of Seurat object}
\item{own_weightage}{ranges from 0 to 1; weight assigned to the cell's own expression in the subtype prediction step}
\item{n_neighbor}{ranges from 5-30; number of neighbors considered to calculate weighted neighbours expression in the subtype prediction step}
}

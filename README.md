# fsctype
Modification of the Original [ScType](https://github.com/IanevskiAleksandr/sc-type) Method for Barcode Annotation using KNN Algorithm


The fsctype method is based on the original `ScType` method used for cluster annotation of single-cell datasets using known markers. 

This variation uses the `data.table` package to substantially speed up calculations and also makes use of the `future` package for parallel execution. 

A typical dataset of 10,000 barcodes or less will run in under 10 seconds in a memory efficient way. You should only rely on parallel execution if your dataset is larger than 50,000 barcodes as the speed up is minimal otherwise and not noticeable. 

This method is different than the original in that instead of relying on precomputed clusters for annotation, it traverses each barcode in a neighborhood graph and uses the k nearest-neighbors to assign a per-cell annotation. 

This method is highly sensitive to the marker genes used. You can use the marker genes from the original ScType database or provide your own, as long as the input marker list is created to look like the output of gene_set_prepare.R from the original method. In other words, you need a nested named list with at least positive markers for your cell types or regions of interest. 

```
pkgs <- c('Seurat', 'dplyr', 'ggplot2', 'data.table')
invisible(lapply(pkgs, require, character.only=TRUE))

source('https://github.com/shahrozeabbas/fsctype/blob/main/R/fsctype.R')

object <- readRDS('seurat_object.rds')
markers <- readRDS('annotation_markers_list.rds')


counts <- as.data.frame(
        object %>% 
            Seurat::ScaleData(verbose=FALSE, features=unique(unlist(markers))) %>% 
            Seurat::GetAssayData(slot='scale.data')
)


cells <- colnames(object)
graph <- object %>% get_graph()
k <- round(sqrt(ncol(object) * 0.5))

predictions <- fsctype(barcodes=cells, graph=graph, counts=counts, markers=markers, n_neighbors=k)


cell_type <- predictions[, prediction]
names(cell_type) <- predictions[, cells]
object <- object %>% AddMetaData(metadata=factor(cell_type), col.name='cell_type') 

g <- object %>% DimPlot(group.by='cell_type', label=TRUE, repel=TRUE, pt.size=1)
ggsave(plot=g, width=12, height=8, filename='plots/pbmc3k_fsctype_umap.png')

```

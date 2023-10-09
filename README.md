# fsctype
### Modification of the Original [ScType](https://github.com/IanevskiAleksandr/sc-type) Method for Barcode Annotation using KNN Algorithm


The fsctype method is based on the original `ScType` method used for cluster annotation of single-cell datasets using known markers. 

This variation uses the `data.table` package to substantially speed up calculations and also makes use of the `future` package for parallel execution. 

A typical dataset of 10,000 barcodes or less will run in under 10 seconds in a memory efficient way. You should only rely on parallel execution if your dataset is larger than 50,000 barcodes as the speed up is minimal otherwise and not noticeable. 

This method is different than the original in that instead of relying on precomputed clusters for annotation, it traverses each barcode in a neighborhood graph and uses the k nearest-neighbors to assign a per-cell annotation. 

This method is highly sensitive to the marker genes used. You can use the marker genes from the original ScType database or provide your own, as long as the input marker list is created to look like the output of `gene_sets_prepare.R` from the original method. In other words, you need a nested named list with at least positive markers for your cell types or regions of interest. 

Your dataset must be processed up to the shared neighborhood graph calculation or some `igraph` object with similarity weights for edges. 

```
pkgs <- c('Seurat', 'dplyr', 'ggplot2', 'data.table')
invisible(lapply(pkgs, require, character.only=TRUE))

object <- CreateSeuratObject(counts=Read10X_h5('filtered_feature_bc_matrix.h5'))

object <- object %>%

                NormalizeData() %>%
                FindVariableFeatures() %>%
                ScaleData(verbose=FALSE) %>%
                RunPCA(verbose=FALSE) %>%
                FindNeighbors(dims=seq(10)) %>%
                RunUMAP(dims=seq(10))
```

Once you have run `FindNeighbors`, you can use this graph downstream for `fsctype` as it relies on the k nearest neighbors of each cell for annotation. The `annotation_markers_list.rds` is a created marker list as defined by `gene_sets_prepare.R`. You can use this from the original `ScType` method or create your own. 

```
source('https://github.com/shahrozeabbas/fsctype/blob/main/R/fsctype.R')

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

ggsave(
  width=12, height=8, filename='fsctype_predictions_umap.png',
  plot=g <- object %>% DimPlot(group.by='cell_type', label=TRUE, repel=TRUE, pt.size=1)
)

```

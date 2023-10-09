require(data.table)

fsctype <- function(barcodes, graph, counts, markers, n_neighbors=20, positive.only=TRUE) {

    scores <- get_scores(counts=counts, markers=markers, positive.only=positive.only)

    neighborhood <- data.table::rbindlist(future.apply::future_sapply(barcodes, function(cell) {
        data.table(cells=get_neighbors(g=graph, cell=cell, k=n_neighbors))
    }, simplify=FALSE), idcol='seed')
    

    predictions <- scores[neighborhood, on='cells', allow.cartesian=TRUE][, 
        .(score=sum(score)), .(seed, type)][, .SD[which.max(score)], by=seed][, 
        .(cells=seed, prediction=type, score)]

    return(predictions)

}

get_graph <- function(object) {
    
    g <- as(object@graphs[[Graphs(object)[2]]], 'dgCMatrix')
    g <- igraph::graph_from_adjacency_matrix(adjmatrix=g, mode='undirected', weighted=TRUE)

    return(g)

}

get_neighbors <- function(g, cell, k=20) {
    
    neighbors <- igraph::ego(g, nodes=igraph::V(g)[cell])[[1]][-1]

    if (length(neighbors) < 1) {
        neighbors <- rep('NA', k)
    } else {
        if (length(neighbors) > k ) {
            neighbors <- neighbors[seq(k)]$name
        } else {
            neighbors <- neighbors$name
        }
    }
    
    return(neighbors)

}

get_scores <- function(counts, markers, positive.only=TRUE) {

    counts <- data.table::melt(
        data.table::as.data.table(copy(counts), keep.rownames='gene_'), 
        id.var='gene_', variable.name='cells', value.name='data'
    )

    marker_stat <- sort(table(unlist(markers)), decreasing=TRUE)

    cell_markers_genes_score <- data.table(
        score_marker_sensitivity=scales::rescale(as.numeric(marker_stat), 
                                    to=c(0, 1), from=c(length(markers), 1)),
        gene_=names(marker_stat), strinmarkersAsFactors=FALSE
    )

    counts[, data := cell_markers_genes_score[counts, on='gene_', score_marker_sensitivity * data]]

    scores <- counts[

        rbindlist(sapply(markers$gs_positive, function(genes) data.table(gene_=genes), simplify=FALSE), idcol='type')
        
        , on='gene_', allow.cartesian=TRUE][, .(score=sum(data) / sqrt(.N)), .(cells, type)]


    if (!positive.only) {

        neg_table <- counts[

        rbindlist(sapply(markers$gs_negative, function(genes) data.table(gene_=genes), simplify=FALSE), idcol='type')
        
        , on='gene_', allow.cartesian=TRUE][, .(neg=sum(data * -1) / sqrt(.N)), .(cells, type)]


        scores <- neg_table[scores, on=c('cells', 'type')][is.na(neg), neg := 0][, .(cells, type, score=score + neg)]
    }
    
    return(scores)

}

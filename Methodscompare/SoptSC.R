library(RSoptSC)
df <- system.file("extdata", "GSE67602_JoostData.csv.bz2", package = "RSoptSC")
gf <- system.file("extdata", "GSE67602_JoostGenes.csv.bz2", package = "RSoptSC")
cf <- system.file("extdata", "GSE67602_JoostCells.csv.bz2", package = "RSoptSC")
af <- system.file("extdata", "GSE67602_JoostAnnotation.csv.bz2", package = "RSoptSC")

GSE67602_Joost <- LoadData(df, gf, cf, af)



logdata <- log10(GSE67602_Joost$data + 1)
gene_expression_threshold <- 0.03
n_features <- 3000
filtered_data <- SelectData(logdata, gene_expression_threshold, n_features)



S <- SimilarityM(lambda = 0.05, 
                 data = filtered_data$M_variable,
                 dims = 3,
                 pre_embed_method = 'tsne',
                 perplexity = 20, 
                 pca_center = TRUE, 
                 pca_scale = TRUE)

low_dim_mapping <- RepresentationMap(similarity_matrix = S$W,
                                     flat_embedding_method = 'tsne',
                                     join_components = FALSE,
                                     perplexity = 35,
                                     theta = 0.5,
                                     normalize = FALSE,
                                     pca = TRUE,
                                     pca_center = TRUE,
                                     pca_scale = TRUE,
                                     dims = 2,
                                     initial_dims = 2)

clusters <- CountClusters(data = S$W, n_comp = 15)
n_clusters <- clusters$upper_bound

plot(c(1:20), 
     clusters$eigs$val[1:20],
     xlab = NA,
     ylab = 'eigenvalues',
     main = 'Eigenvalues of the Graph Laplacian')



output_NMF <- NMF::nmf(x = Matrix::as.matrix(S$W),
                       rank = n_clusters,
                       method = 'lee',
                       seed = 'nndsvd',
                       .options = 'nP');
H <- NMF::basis(output_NMF)


labels <- apply(H, 1, function(x){
  which(x == max(x))})



# define a scheme for the coloring.  This is a required parameter for plotting discrete (factor) data
colorscale <- ColorHue(n = length(unique(labels)))
colorscale <- colorscale$hex.1.n.

# plot clusters
FeatureScatterPlot(flat_embedding = low_dim_mapping$flat_embedding,
                   feature = as.factor(labels),
                   title = "NMF Cluster Labeling",
                   subtitle = "t-SNE Embedding",
                   featurename = "Cluster ID",
                   colorscale = colorscale)


true_labels <- GSE67602_Joost$annotation
# plot clusters
FeatureScatterPlot(flat_embedding = low_dim_mapping$flat_embedding,
                   feature = as.factor(true_labels),
                   title = "True Labeling",
                   subtitle = "t-SNE Embedding",
                   featurename = "Annotated Cell Types")


markers <- GetMarkerTable(counts_data = filtered_data,
                          cluster_labels = labels,
                          H = H,
                          n_sorted = 25,
                          
                          gene_expression_threshold = 3,

                          use_H = TRUE)

PlotTopN(GSE67602_Joost$data = log(GSE67602_Joost$data+1),
         gene_names = gene_names,
         cluster_labels = labels,
         markers = markers$all,
         n_features = 5)


cluster_ptime <- FindRootCluster(cluster_labels = labels,
                                 flat_embedding = low_dim_mapping$flat_embedding,
                                 dist_graph = low_dim_mapping$dist_graph,
                                 dist_flat = low_dim_mapping$dist_flat,
                                 reverse = TRUE)


root_cell <- FindRootCell(use_flat_dist = FALSE,
                          cluster_order_by = "distance",
                          cell_order_by = "distance",
                          graph_cluster_mst = cluster_ptime$cluster_mst,
                          dist_graph  = low_dim_mapping$dist_graph,
                          dist_flat = low_dim_mapping$dist_flat,
                          cluster_labels = labels,
                          root_cluster = cluster_ptime$root_cluster)


cluster_predecessors <- GetPredecessors(cluster_ptime$cluster_mst, cluster_ptime$root_cluster)
cluster_dtree <- GetDominatorTree(cluster_predecessors, cluster_ptime$graph_cluster)
PlotLineage(cluster_dtree)

pseudotime <- low_dim_mapping$dist_graph[root_cell,]

# plot pseudotime
FeatureScatterPlot(flat_embedding = low_dim_mapping$flat_embedding,
                   feature = pseudotime,
                   title = "Pseudotime Labeling",
                   subtitle = "t-SNE Embedding",
                   featurename = "Pseudotime Distance")


gene_index <- which(gene_names == 'Krt14')
data_gene <- logdata[gene_index,]

# plot features
FeatureScatterPlot(flat_embedding = low_dim_mapping$flat_embedding,
                   feature = data_gene,
                   title = "Krt14 Expression",
                   subtitle = "t-SNE Embedding",
                   featurename = "Log Expression")


gene_index <- which(gene_names == 'Krt10')
data_gene <- logdata[gene_index,]

# plot features
FeatureScatterPlot(flat_embedding = low_dim_mapping$flat_embedding,
                   feature = data_gene,
                   title = "Krt10 Expression",
                   subtitle = "t-SNE Embedding",
                   featurename = "Log Expression")


gene_index <- which(gene_names == 'Lor')
data_gene <- logdata[gene_index,]

# plot features
FeatureScatterPlot(flat_embedding = low_dim_mapping$flat_embedding,
                   feature = data_gene,
                   title = "Lor Expression",
                   subtitle = "t-SNE Embedding",
                   featurename = "Log Expression")

plot(ViolinPlotExpression(data = logdata,
                          gene_names = gene_names,
                          labels = labels,
                          gene_name = "Krt14"))

plot(ViolinPlotExpression(data = logdata,
                          gene_names = gene_names,
                          labels = labels,
                          gene_name = "Krt10"))

plot(ViolinPlotExpression(data = logdata,
                          gene_names = gene_names,
                          labels = labels,
                          gene_name = "Lor"))

library(knitr)
library(kableExtra)

lig_rec_path <- system.file("extdata", "tgfb_lig_rec.tsv", package = "RSoptSC")
rec_target_path <- system.file("extdata", "tgfb_rec_target_both.tsv", package = "RSoptSC")
pathway <- ImportPathway(lig_table_path = lig_rec_path,
                         rec_table_path = rec_target_path,
                         data = logdata,
                         gene_names = gene_names)



pathway$pathway %>% kable() %>% kable_styling()

Pmats <- GetSignalingPartners(logdata,
                              gene_names,
                              pathway$pathway_removed)

SigPlot(P = Pmats$P_agg,
        cluster_labels = labels,
        lig_cells_per_cluster = 5,
        rec_cells_per_cluster = 5,
        rgb_gap = 0.06,
        title_text = "Signaling Between Cells")


cluster_P <- ClusterSig(Pmats$P_agg,
                        cluster_labels = labels)
SigPlot(P = cluster_P,
        cluster_labels = c(1:n_clusters),
        rgb_gap = 0.60,
        title_text = "Signaling Between Clusters")

gene_list = c('Tgfb1', 'Tgfb2', 
              'Tgfbr1', 'Tgfbr2',
              'Zeb2','Smad2','Wnt4','Wnt11','Bmp7','Sox9','Notch1')

PlotClusterExpression(data = logdata,
                      gene_names = gene_names,
                      cluster_labels = labels, 
                      markers = gene_list)


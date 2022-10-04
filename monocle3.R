
#BiocManager::install("monocle")#,force = TRUE)
library(Seurat)
library(monocle3)
#devtools::install_github('cole-trapnell-lab/monocle3')

#library(VGAM)
setwd("/home/oulas/scRNA-Seq/")
#library(cellrangerRkit)


#dat <- readRDS("COVID/COVID-19_brain_snRNA-seq_choroid_plexus_final_seurat_v3.2.3.rds")
#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(parenchyma_cortex@assays$RNA@data), 'sparseMatrix')
#data <- as(as.matrix(LAM1@assays$RNA@data), 'sparseMatrix')

pData<-parenchyma_cortex@meta.data
#pd <- new('AnnotatedDataFrame', data = parenchyma_cortex@meta.data)
#pd <- new('AnnotatedDataFrame', data = LAM1@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
#fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
monocle_cds <- new_cell_data_set(expression_data = data,cell_metadata = pData,gene_metadata =fData)#,cell_metadata = pd, gene_metadata = fd)



pData(monocle_cds)
fData(monocle_cds)

#Normalize and Pre-process the data
monocle_cds <- preprocess_cds(monocle_cds, num_dim = 100)

## Step 2: Remove batch effects with cell alignment
monocle_cds <- align_cds(monocle_cds, alignment_group = "Biogroup")

## Step 3: Reduce the dimensions using UMAP
monocle_cds <- reduce_dimension(monocle_cds, reduction_method="UMAP")
plot_cells(monocle_cds,color_cells_by="cellID")


## Step 4: Cluster the cells
monocle_cds <- cluster_cells(monocle_cds, resolution=1e-5)
plot_cells(monocle_cds,color_cells_by="Biogroup",label_principal_points = T)


monocle_cds_subset <- choose_cells(monocle_cds)
  



#Order cells in pseudotime along a trajectory
## Step 5: Learn a graph first then subset to see trajectories
monocle_cds <- learn_graph(monocle_cds)
#monocle_cds_subset<- learn_graph(monocle_cds_subset)

#subsetting
monocle_cds_subset <- monocle_cds[,pData(monocle_cds)$cellID2 == "Astrocyte"]
monocle_cds_subset <- reduce_dimension(monocle_cds_subset, reduction_method="UMAP")
plot_cells(monocle_cds_subset,color_cells_by="Biogroup")
monocle_cds_subset <- cluster_cells(monocle_cds_subset, resolution=1e-5)
plot_cells(monocle_cds_subset,color_cells_by="Biogroup")

plot_cells(monocle_cds_subset,color_cells_by="Biogroup",label_principal_points = T)

# a helper function to identify the root principal points: does not work
get_correct_root_state <- function(cds, cell_phenotype, root_type){
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
  
  closest_vertex <-
    cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    V(cds@minSpanningTree)$name[as.numeric(names
                                           (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin){
  cell_ids <- which(colData(cds)[, "cellID"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}


root_node_ids = get_correct_root_state(monocle_cds,"cellID", "Astrocyte")
cellids<-unique(pData(monocle_cds)$cellID)
root_pr_nodes<-c()
for(i in 1:length(cellids)){
  root_pr_nodes=c(root_pr_nodes,get_earliest_principal_node(monocle_cds,cellids[i]))
}

#plot_cell_trajectory(monocle_cds_subset)
monocle_cds_subset <- order_cells(monocle_cds_subset,root_pr_nodes = c("Y_326"))#c("Y_179","Y_87","Y_15")

## Step 6: Order cells
monocle_cds <- order_cells(monocle_cds,root_pr_nodes=root_pr_nodes)#,root_pr_nodes =  "Astrocyte")

plot_cells(monocle_cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)


plot_cells(monocle_cds_subset,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)




#Finding genes that change as a function of pseudotime
#Identifying the genes that change as cells progress along a trajectory is a core objective of this 
#type of analysis

#We can use the graph_test(), passing it neighbor_graph="principal_graph", which tells it to test 
#whether cells at similar positions on the trajectory have correlated expression
monocle_cds_subset_pr_test_res <- graph_test(monocle_cds_subset, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(monocle_cds_subset_pr_test_res[order(monocle_cds_subset_pr_test_res$q_value),], q_value < 0.05))

#plot interesting genes that score as highly significant according to graph_test()
plot_cells(monocle_cds_subset, genes=pr_deg_ids[4:8],
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)






#we can collect the trajectory-variable genes into modules:
gene_module_df <- find_gene_modules(monocle_cds_subset[pr_deg_ids,], resolution=c(10^seq(-6,-1)))



#Here we plot the aggregate module scores within each biogroup Case vs. Control
cell_group_df <- tibble::tibble(cell=row.names(colData(monocle_cds_subset)), 
                                  cell_group=colData(monocle_cds_subset)$Biogroup)#
agg_mat <- aggregate_gene_expression(monocle_cds_subset, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")






#We can also pass gene_module_df to plot_cells() .

plot_cells(monocle_cds_subset,
           genes=gene_module_df %>% filter(module %in% c(20, 10, 7, 1)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)



#Monocle offers another plotting function that can sometimes give a clearer view of a gene's dynamics 
#along a single path. You can select a path with choose_cells() or by subsetting the cell data set by cluster, 
#cell type, or other annotation that's restricted to the path. Let's pick one such path for some Astro cells

selected_genes <- pr_deg_ids[13:18]
selected_genes_lineage_monocle_cds_subset <- monocle_cds_subset[rowData(monocle_cds_subset)$gene_short_name %in% selected_genes,
                       colData(monocle_cds_subset)$cellID2 %in% c("Astrocyte")]


#The function plot_genes_in_pseudotime() takes a small set of genes and shows you 
#their dynamics as a function of pseudotime:
#You can see that some genes may be activated before others.
plot_genes_in_pseudotime(selected_genes_lineage_monocle_cds_subset,
                         color_cells_by="Biogroup",
                         min_expr=0.5)







colnames(pData(monocle_cds))

pData(monocle_cds)$Total_mRNAs <- Matrix::colSums(exprs(monocle_cds))

monocle_cds <- monocle_cds[,pData(monocle_cds)$Total_mRNAs < 1e6]

upper_bound <- 10^(mean(log10(pData(monocle_cds)$Total_mRNAs)) +
                     2*sd(log10(pData(monocle_cds)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(monocle_cds)$Total_mRNAs)) -
                     2*sd(log10(pData(monocle_cds)$Total_mRNAs)))

qplot(Total_mRNAs, data = pData(monocle_cds),  geom =
        "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)


L <- log(exprs(monocle_cds[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
melted_dens_df <- reshape2::melt(Matrix::t(scale(Matrix::t(L)))) #takes too long

# Plot the distribution of the standardized gene expression values. takes too long
qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")


#Classifying cells by type Recommended
TAGLN_id <- row.names(subset(fData(monocle_cds), gene_short_name == "TAGLN"))
NRF2_id <- row.names(subset(fData(monocle_cds), gene_short_name == "NFE2L2"))
ACTA2_id <- row.names(subset(fData(monocle_cds),
                             gene_short_name == "ACTA2"))

cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "SmoothMC", classify_func =
                     function(x) { x[NRF2_id,] >= 1 & x[ACTA2_id,] > 1 })

cth <- addCellType(cth, "SmoothMC", classify_func =
                     function(x) { x[NRF2_id,] >= 1 })
#cth <- addCellType(cth, "SmoothMC2", classify_func = function(x)
#{ x[TAGLN_id,] < 1 & x[ACTA2_id,] > 1 })

monocle_cds <- classifyCells(monocle_cds, cth, 0.1)
table(pData(monocle_cds)$CellType)

pie <- ggplot(pData(monocle_cds),
              aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())


#Clustering cells without marker genes Alternative
#Monocle provides an algorithm you can use to impute the types of the "Unknown" cells
disp_table <- dispersionTable(monocle_cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
monocle_cds <- setOrderingFilter(monocle_cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(monocle_cds)

#Now we're ready to try clustering the cells:
# monocle_cds@auxClusteringData[["tSNE"]]$variance_explained <- NULL
#plot_pc_variance_explained(monocle_cds, return_all = F) # norm_method='log'

monocle_cds <- reduceDimension(monocle_cds, max_components = 2, num_dim = 6,
                               reduction_method = 'tSNE', verbose = T)

monocle_cds <- clusterCells(monocle_cds, num_clusters = 8)
plot_cell_clusters(monocle_cds, 1, 2,color = "Disease",
                   markers = c("NFE2L2"))


plot_cell_clusters(monocle_cds, 1, 2,color = "cellID")
plot_cell_clusters(monocle_cds, 1, 2,color = "Cluster")

#Clustering cells using marker genes Recommendedcd-hit.R
#First, we'll select a different set of genes to use for clustering the cells. Before we just picked genes that 
#were highly expressed and highly variable. Now, we'll pick genes that co-vary with our markers. In a sense,
#we'll be building a large list of genes to use as markers, so that even if a cell doesn't have TAGLN, 
#it might be recognizable as a Smooth muscle based on other genes.
# marker_diff <- markerDiffTable(monocle_cds[expressed_genes,],
#                                cth,
#                                #residualModelFormulaStr = "~num_genes_expressed",
#                                cores = 1)

#Monocle provides a handy function for ranking genes by how restricted their expression is for each type.
# candidate_clustering_genes <-c("ACTA2","TAGLN")
# row.names(subset(marker_diff, qval < 0.01))
# marker_spec <-
#   calculateMarkerSpecificity(monocle_cds[candidate_clustering_genes,], cth)
# head(selectTopMarkers(marker_spec, 3))
# 
# #To cluster the cells, we'll choose the top 500 markers for each of these cell types:
# semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 2)$gene_id)#500
# monocle_cds <- setOrderingFilter(monocle_cds, semisup_clustering_genes)
# plot_ordering_genes(monocle_cds)
# monocle_cds <- reduceDimension(monocle_cds, max_components = 2, num_dim = 3,
#                                norm_method = 'log',
#                                reduction_method = 'tSNE',
#                                #residualModelFormulaStr = "~CellType + num_genes_expressed",
#                                verbose = T)
# monocle_cds <- clusterCells(monocle_cds, num_clusters = 2)
# plot_cell_clusters(monocle_cds, 1, 2, color = "Disease")


monocle_cds_subset <- monocle_cds[,pData(monocle_cds)$cellID2 == "Astrocyte"]
#monocle_cds_subset <- monocle_cds[,pData(monocle_cds)$Cluster == "6"]

#Trajectory step 1: choose genes that define a cell's progress
#One effective way to isolate a set of ordering genes is to 
#simply compare the cells collected at the beginning of the process to 
#those at the end and find the differentially expressed genes, as described above. 
#The command below will find all genes that are differentially expressed in response to the switch from control to disease:
diff_test_res <- differentialGeneTest(monocle_cds_subset[expressed_genes,],
                                      fullModelFormulaStr = "~Biogroup")

diff_test_res_sorted <- diff_test_res[order(diff_test_res$qval,decreasing = F),]

plot_genes_jitter(monocle_cds_subset[rownames(diff_test_res_sorted)[1:3],],
                  grouping = "Biogroup",
                  min_expr = 0.1)

cds_subsetTopGenes <- monocle_cds_subset[rownames(diff_test_res_sorted)[1:3],]
plot_genes_in_pseudotime(cds_subsetTopGenes, color_by = "Biogroup")

ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

#Once we have a list of gene ids to be used for ordering, we need to set them in the monocle_cds_subset object, 
#because the next several functions will depend on them.

monocle_cds_subset <- setOrderingFilter(monocle_cds_subset, ordering_genes)
plot_ordering_genes(monocle_cds_subset)

#Trajectory step 2: reduce data dimensionality
#Next, we will reduce the space down to one with two dimensions, which we will be able to easily visualize and interpret while Monocle is ordering the cells.

#Not sure if this is needed...
monocle_cds_subset <- estimateDispersions(monocle_cds_subset)
monocle_cds_subset <- reduceDimension(monocle_cds_subset, max_components=2, method = 'DDRTree',ncenter = 50,auto_param_selection = F)


#Trajectory step 3: order cells along the trajectory
#Now that the space is reduced, it's time to order the cells using the orderCells function as shown below
monocle_cds_subset <- orderCells(monocle_cds_subset)




#The trajectory has a tree-like structure. We can see that the cells collected at time zero are located near one of 
#the tips of the tree, while the others are distributed amongst the two "branches". Monocle doesn't know a priori which 
#of the trajectory of the tree to call the "beginning", so we often have to call orderCells again using the root_state argument to specify the beginning. 
#First, we plot the trajectory, this time coloring the cells by "State":
plot_cell_trajectory(monocle_cds_subset, color_by = "State")
plot_cell_trajectory(monocle_cds_subset, color_by = "Biogroup")
plot_cell_trajectory(monocle_cds_subset, color_by = "Disease")
# plot_cell_trajectory(monocle_cds_subset, color_by = "State") +
#   facet_wrap(~State, nrow = 1)

plot_cell_trajectory(monocle_cds_subset, color_by = "Pseudotime")

#"State" is just Monocle's term for the segment of the tree. The function below is handy 
#for identifying the State which contains most of the cells from Biogroup Control. We can then pass that to orderCells:

GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$Biogroup)[,"Control"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
#Can view object
monocle_cds_subset$Pseudotime

monocle_cds_subset <- orderCells(monocle_cds_subset, root_state = GM_state(monocle_cds_subset))
plot_cell_trajectory(monocle_cds_subset, color_by = "Pseudotime")

#Finding Genes that Change as a Function of Pseudotime
diff_test_resPseudo <- differentialGeneTest(cds_subsetTopGenes,
                                            fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_test_resPseudo[,c("gene_short_name", "pval", "qval")]

plot_genes_in_pseudotime(cds_subsetTopGenes, color_by = "Disease")


#Clustering Genes by Pseudotemporal Expression Pattern
plot_pseudotime_heatmap(monocle_cds_subset[rownames(diff_test_resPseudo),],
                        num_clusters = 2,
                        cores = 1,
                        show_rownames = T)

#Analyzing Branches in Single-Cell Trajectories
#BEAM takes as input a CellDataSet that's been ordered with orderCells and the name of a branch point in the trajectory. It returns a table of 
#significance scores for each gene. Genes that score significant are said to be branch-dependent in their expression.

BEAM_res <- BEAM(monocle_cds_subset, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

plot_genes_branched_heatmap(monocle_cds_subset[row.names(subset(BEAM_res,
                                                                qval < 1e-4)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)

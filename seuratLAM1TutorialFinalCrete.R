library(dplyr)
library(Seurat)
library(patchwork)
library(multtest)
library(metap)
library(ggplot2)
library(cowplot)
library(enrichR)

setwd("/home/oulas/scRNA-Seq")

# Load the dataset
LAM1.data <- Read10X(data.dir = "../GSE135851_RAW/LAM1/")

# Initialize the Seurat object with the raw (non-normalized data).
LAM1 <- CreateSeuratObject(counts = LAM1.data, project = "LAM1", min.cells = 3, min.features = 200)
LAM1


#QC and selecting cells for further analysis
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
LAM1[["percent.mt"]] <- PercentageFeatureSet(LAM1, pattern = "^MT-")
head(LAM1@meta.data, 5)


# Visualize QC metrics as a violin plot
VlnPlot(LAM1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(LAM1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(LAM1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Filtering data
LAM1 <- subset(LAM1, subset = nFeature_RNA > 100 & nFeature_RNA < 6500 & percent.mt < 75)


#Normalizing the data. Normalized values are stored in LAM1[["RNA"]]@data. *
LAM1 <- NormalizeData(LAM1, normalization.method = "LogNormalize", scale.factor = 10000)


#Identification of highly variable features (feature selection) *
LAM1 <- FindVariableFeatures(LAM1, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(LAM1), 10)


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(LAM1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#Scaling the data

all.genes <- rownames(LAM1)
LAM1 <- ScaleData(LAM1, features = all.genes)


#Perform linear dimensional reduction
LAM1 <- RunPCA(LAM1, features = VariableFeatures(object = LAM1))


# Examine and visualize PCA results a few different ways
print(LAM1[["pca"]], dims = 1:5, nfeatures = 5)


VizDimLoadings(LAM1, dims = 1:2, reduction = "pca")


DimPlot(LAM1, reduction = "pca")


DimHeatmap(LAM1, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(LAM1, dims = 1:15, cells = 500, balanced = TRUE)


#Determine the dimensionality of the dataset

# NOTE: JackStraw process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
#LAM1 <- JackStraw(LAM1, num.replicate = 100) #approx 3 mins
#LAM1 <- ScoreJackStraw(LAM1, dims = 1:20)
#JackStrawPlot(LAM1, dims = 1:15)
ElbowPlot(LAM1)

#Cluster the cells

LAM1 <- FindNeighbors(LAM1, dims = 1:10)
LAM1 <- FindClusters(LAM1, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(LAM1), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# umap-learn)
LAM1 <- RunUMAP(LAM1, dims = 1:10)
LAM1 <- RunTSNE(LAM1,dims = 1:10)

# note that you can set label = TRUE or use the LabelClusters function to help label
# individual clusters
DimPlot(LAM1, reduction = "umap")
DimPlot(LAM1, reduction = "tsne")


#Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 2
cluster2.markers <- FindMarkers(LAM1, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(LAM1, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones - approx 4.5 mins with default test - approx 12 mins with MAST
#BiocManager::install("MAST")

LAM1.markers <- FindAllMarkers(LAM1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)#,test.use = "MAST"

saveRDS(LAM1.markers, file = "LAM1.markers.rds")
#LAM1.markers<-readRDS(file = "LAM1.markers.rds")

LAM1.markerstop2<-LAM1.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

LAM1.markerstop100<-LAM1.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC)


#Visualize markers
VlnPlot(LAM1, features = LAM1.markerstop2$gene[LAM1.markerstop2$cluster==16])


FeaturePlot(LAM1, features = LAM1.markerstop2$gene[LAM1.markerstop2$cluster==16])


DoHeatmap(LAM1, features = LAM1.markerstop2$gene) + NoLegend()


#Assigning cell type identity to clusters
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

dbstouse <- c("PanglaoDB_Augmented_2021","CellMarker_Augmented_2021")
annotatedclusters<-c()
listofresults<-c()
for(m in as.numeric(as.character(unique(LAM1.markerstop100$cluster)))){
  markers_to_search<-LAM1.markerstop100$gene[LAM1.markerstop100$cluster==m]
  if (websiteLive) {
    enriched <- enrichr(markers_to_search, dbstouse)
  }
  celltype<-c(enriched[["PanglaoDB_Augmented_2021"]]$Term[1],enriched[["CellMarker_Augmented_2021"]]$Term[1])
  print(paste("Cluster",m,celltype,sep = " "))
  names(enriched)<-paste(names(enriched),m,sep = "_")
  listofresults<-c(listofresults,enriched)
  annotatedclusters<-c(annotatedclusters,toString(celltype[1]))
}
if (websiteLive) plotEnrich(listofresults[["PanglaoDB_Augmented_2021_16"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "PanglaoDB")+plotEnrich(listofresults[["CellMarker_Augmented_2021_16"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CellMarker")
print(annotatedclusters)

#Change labels
new.cluster.ids <- make.names(annotatedclusters,unique = T)
names(new.cluster.ids) <- levels(LAM1)
LAM1 <- RenameIdents(LAM1, new.cluster.ids)
unique(Idents(LAM1))
DimPlot(LAM1, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
##############################################################################################################################################


##########################################Run Integrated analysis#############################################################################
LAM1.data <- Read10X(data.dir = "../GSE135851_RAW/LAM1/")
WT.data <- Read10X(data.dir = "../GSE135851_RAW/Donor1/")

# Initialize the Seurat object with the raw (non-normalized data).
LAM1 <- CreateSeuratObject(counts = LAM1.data, project = "LAM1", min.cells = 3, min.features = 200)
LAM1
WT <- CreateSeuratObject(counts = WT.data, project = "WT", min.cells = 3, min.features = 200)
WT

LAM1[["percent.mt"]] <- PercentageFeatureSet(LAM1, pattern = "^MT-")
WT[["percent.mt"]] <- PercentageFeatureSet(WT, pattern = "^MT-")
LAM1_WT.list<-c(LAM1,WT)


LAM1_WT.list <- lapply(X = LAM1_WT.list, FUN = function(x) {
  #Filtering data
  #x <- subset(x, subset = nFeature_RNA >= 0 & nFeature_RNA < 9500 & percent.mt < 100)
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


features <- SelectIntegrationFeatures(object.list = LAM1_WT.list)
LAM1_WT.anchors <- FindIntegrationAnchors(object.list = LAM1_WT.list, dims = 1:20)#anchor.features = features
LAM1_WT.combined <- IntegrateData(anchorset = LAM1_WT.anchors, dims = 1:20)
#saveRDS(LAM1_WT.combined, file = "LAM1_WT.combined.rds")
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(LAM1_WT.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
LAM1_WT.combined <- ScaleData(LAM1_WT.combined, verbose = FALSE)
LAM1_WT.combined <- RunPCA(LAM1_WT.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
LAM1_WT.combined <- RunUMAP(LAM1_WT.combined, reduction = "pca", dims = 1:20)
LAM1_WT.combined <- RunTSNE(LAM1_WT.combined,reduction = "pca",dims = 1:10)
#Louvein method or improved Leiden  ??
LAM1_WT.combined <- FindNeighbors(LAM1_WT.combined, reduction = "pca", dims = 1:20)
LAM1_WT.combined <- FindClusters(LAM1_WT.combined, resolution  = 0.3)
#Can view object
LAM1_WT.combined[[]]
#LAM1_WT.combined[['groups']] <- sample(x = c('g1', 'g2'), size = ncol(x = LAM1_WT.combined), replace = TRUE)
p1combined <- DimPlot(LAM1_WT.combined, reduction = "umap", group.by = "orig.ident")
p2combined <- DimPlot(LAM1_WT.combined, reduction = "umap", label = TRUE)
p1combined+p2combined



DimPlot(LAM1_WT.combined, reduction = "umap", split.by = "orig.ident")
DimPlot(LAM1_WT.combined, reduction = "tsne", split.by = "orig.ident")


#Identify conserved cell type markers
# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(LAM1_WT.combined) <- "RNA"



LAM1_WT.markers <- FindAllMarkers(LAM1_WT.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(LAM1_WT.markers, file = "LAM1_WT.markers.rds")
LAM1_WT.markerstop1<-LAM1_WT.markers %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC) #100 gives best results

LAM1_WT.markerstop2<-LAM1_WT.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC) #100 gives best results

LAM1_WT.markerstop100<-LAM1_WT.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC) #100  gives best results

LAM1_WT.markerstop200<-LAM1_WT.markers %>%
  group_by(cluster) %>%
  slice_max(n = 200, order_by = avg_log2FC) #200 gives second best results


FeaturePlot(LAM1_WT.combined, features = LAM1_WT.markerstop1$gene[c(1,2,3,4,5,6,9,10,21)])
FeaturePlot(LAM1_WT.combined, features = LAM1_WT.markers$gene[c(6930)])


#Run Enrich R on to 100 Markers
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

annotatedclusters<-c()
listofresults<-c()
dbstouse <- c("PanglaoDB_Augmented_2021","CellMarker_Augmented_2021")
for(m in as.numeric(as.character(unique(LAM1_WT.markerstop100$cluster)))){
  markers_to_search<-LAM1_WT.markerstop100$gene[LAM1_WT.markerstop100$cluster==m]
  if (websiteLive) {
    enriched <- enrichr(markers_to_search, dbstouse)
  }
  celltype<-c(enriched[[1]]$Term[1],enriched[[2]]$Term[1])
  print(paste("Cluster",m,celltype,sep = " "))
  annotatedclusters<-c(annotatedclusters,toString(celltype[1]))
  names(enriched)<-paste(names(enriched),m,sep = "_")
  listofresults<-c(listofresults,enriched)
}

if (websiteLive) plotEnrich(listofresults[["PanglaoDB_Augmented_2021_12"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "PanglaoDB")+plotEnrich(listofresults[["CellMarker_Augmented_2021_12"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CellMarker")

new.cluster.ids <- make.names(annotatedclusters,unique = T)
names(new.cluster.ids) <- levels(LAM1_WT.combined)
LAM1_WT.combined <- RenameIdents(LAM1_WT.combined, new.cluster.ids)


DimPlot(LAM1_WT.combined, reduction = "umap",label = TRUE)

markers.to.plot <- unique(LAM1_WT.markerstop2$gene) 
DotPlot(LAM1_WT.combined, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "orig.ident") +
  RotatedAxis()
DotPlot(LAM1_WT.combined, features = "NFE2L2", cols = c("blue", "red"), dot.scale = 8, split.by = "orig.ident") +
  RotatedAxis()


#Change the Labels
LAM1_WT.combined$celltype.orig.ident <- paste(Idents(LAM1_WT.combined), LAM1_WT.combined$orig.ident, sep = "_")
LAM1_WT.combined$celltype <- Idents(LAM1_WT.combined)
Idents(LAM1_WT.combined) <- "celltype.orig.ident"
unique(Idents(LAM1_WT.combined))

LAM_WT.markersSM <- FindMarkers(LAM1_WT.combined, ident.1 = "Smooth.Muscle.Cells_LAM1", ident.2 = "Smooth.Muscle.Cells_WT", verbose = FALSE)
LAM_WT.markersSM <- FindMarkers(LAM1_WT.combined, ident.1 = "Endothelial.Cells.1_LAM1", ident.2 = "Endothelial.Cells.1_WT", verbose = FALSE)
indexsigclSM<-which(LAM_WT.markersSM$p_val_adj<=0.05)
LAM_WT.markersSM_sig<-LAM_WT.markersSM[indexsigclSM,]

Up<-rownames(LAM_WT.markersSM_sig[which(LAM_WT.markersSM_sig$avg_log2FC>=0),])
Down<-rownames(LAM_WT.markersSM_sig[which(LAM_WT.markersSM_sig$avg_log2FC<0),])

dbstouse <- c("Old_CMAP_up","Old_CMAP_down","LINCS_L1000_Chem_Pert_up","LINCS_L1000_Chem_Pert_down","LINCS_L1000_Ligand_Perturbations_down","LINCS_L1000_Ligand_Perturbations_up") #PanglaoDB_Augmented_2021
dbstousePaths<-c("KEGG_2021_Human","MSigDB_Hallmark_2020","GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015","ProteomicsDB_2020","PanglaoDB_Augmented_2021","CellMarker_Augmented_2021","Human_Gene_Atlas")

if(length(Up) > 200){
  Up<-Up[1:200]
}
if(length(Down) > 200){
  Down<-Down[1:200]
}

if (websiteLive) {
  enrichedUp <- enrichr(Up, dbstouse)
}

if (websiteLive) {
  enrichedDown <- enrichr(Down, dbstouse)
}
#names(enriched)
#celltype<-c(enriched[[1]]$Term[1],enriched[[2]]$Term[1])
#if (websiteLive) plotEnrich(enriched[["Old_CMAP_up"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CMAP ASTRO UP")+plotEnrich(enriched[["Old_CMAP_down"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CMAP ASTRO Down")
#if (websiteLive) plotEnrich(enriched[["KEGG_2021_Human"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "KEGG")+plotEnrich(enriched[["MSigDB_Hallmark_2020"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "MSigDB")+plotEnrich(enriched[["GO_Molecular_Function_2015"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "GO Molecular Function")

if (websiteLive) {
  enrichedPaths <- enrichr(c(Up,Down), dbstousePaths)
}
if (websiteLive) plotEnrich(enrichedPaths[["KEGG_2021_Human"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = paste("KEGG","LAM SM",sep=" "))#+plotEnrich(enriched[["MSigDB_Hallmark_2020"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "MSigDB")+plotEnrich(enriched[["GO_Molecular_Function_2015"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "GO Molecular Function")

pathsSM<-enrichedPaths[["KEGG_2021_Human"]]

drugsUp<-enrichedUp[["Old_CMAP_down"]]$Term[which(enrichedUp[["Old_CMAP_down"]]$P.value <0.05)]
drugsUp<-gsub("-[0-9]+$","",drugsUp)

drugsDown<-enrichedDown[["Old_CMAP_up"]]$Term[which(enrichedDown[["Old_CMAP_up"]]$P.value <0.05)]
drugsDown<-gsub("-[0-9]+$","",drugsDown)
if (websiteLive) plotEnrich(enrichedUp[["Old_CMAP_down"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CMAP SM UP")+plotEnrich(enrichedDown[["Old_CMAP_up"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CMAP SM Down")

drugInfo<-read.delim("COVID/drug_repurposing_hub.txt")
AlldrugsCellIDs<-c(drugsUp,drugsDown)
foundds<-0
MOA<-c()
for(d in 1:length(AlldrugsCellIDs)){
  indexfoundd<-grep(AlldrugsCellIDs[d],drugInfo$pert_iname)
  if(!length(indexfoundd)==0){
    foundds<-foundds+1
    #MOA<-c(MOA,drugInfo$moa[indexfoundd])
    MOA<-rbind(MOA,cbind(drugInfo$pert_iname[indexfoundd],drugInfo$moa[indexfoundd],drugInfo$disease_area[indexfoundd]))
  }
}
indexesnull1<-which(MOA[,2] == "")
if(!length(indexesnull1)==0){
  MOA<-MOA[-indexesnull1,]
}
indexesnull2<-which(MOA[,3] == "")
if(!length(indexesnull2)==0){
  MOA<-MOA[-indexesnull2,]
}

# sort by Freq
tableMOAs<-as.data.frame(table(MOA[,3]))
tableMOAs <- tableMOAs[order(tableMOAs$Freq,decreasing = T),]
write.table(tableMOAs,"tableMOAs_LAM_SM.txt",quote = F,row.names = F,sep = "\t")

# if (websiteLive) {
#   enriched <- enrichr(Up[1:200], dbstouse)
# }
# names(enriched)
# celltype<-c(enriched[[1]]$Term[1],enriched[[2]]$Term[1])
# if (websiteLive) plotEnrich(enriched[["Old_CMAP_up"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CMAP UP")+plotEnrich(enriched[["Old_CMAP_down"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CMAP Down")
# if (websiteLive) plotEnrich(enriched[["KEGG_2021_Human"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "KEGG")+plotEnrich(enriched[["MSigDB_Hallmark_2020"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "MSigDB")+plotEnrich(enriched[["GO_Molecular_Function_2015"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "GO Molecular Function")

#Subsetting cluster
#Identify differential expressed genes across conditions
LAM1_WT.combinedclSM<-subset(LAM1_WT.combined, idents = c("Smooth.Muscle.Cells_LAM1","Smooth.Muscle.Cells_WT"))
LAM1_WT.combinedclSM <- subset(LAM1_WT.combined, idents = c("Endothelial.Cells.1_LAM1","Endothelial.Cells.1_WT"))
#Idents(LAM1_WT.combinedcl9)<-"LAM1"
LAM1_WT.combinedclSM<-ScaleData(LAM1_WT.combinedclSM, features = rownames(LAM_WT.markersSM_sig))

theme_set(theme_cowplot())
avg.clSM <- as.data.frame(log1p(AverageExpression(LAM1_WT.combinedclSM, verbose = FALSE)$RNA))
avg.clSM$gene <- rownames(avg.clSM)
genes.to.label = rownames(LAM_WT.markersSM_sig)[1:10]
avg.clSM$gene[which(!avg.clSM$gene%in%genes.to.label)]<-""

p1 <- ggplot(avg.clSM, aes(`Smooth.Muscle.Cells_LAM1`, `Smooth.Muscle.Cells_WT`,label=gene)) + geom_point() + ggtitle("LAM1 signature genes")+geom_text()
p1

DoHeatmap(LAM1_WT.combinedclSM, features = c(Up[1:50],Down),raster=F,size = 4)+
  theme(axis.text.y = element_text(size = 8))

DoHeatmap(LAM1_WT.combined, features = "NFE2L2",raster=F,size = 4)+
  theme(axis.text.y = element_text(size = 8))

#Simulation Bulk RNA-seq
LAM1Cells<-unique(Idents(LAM1_WT.combined))[grep("LAM1",unique(Idents(LAM1_WT.combined)))]
WTCells<-unique(Idents(LAM1_WT.combined))[grep("WT",unique(Idents(LAM1_WT.combined)))]
LAM_WT.markersbulk <- FindMarkers(LAM1_WT.combined, ident.1 = LAM1Cells,ident.2 = WTCells, verbose = FALSE)
                                    #c('Dendritic.Cells_LAM1','Dendritic.Cells.1_LAM1','Pulmonary.Alveolar.Type.II.Cells.2_LAM1','Pulmonary.Alveolar.Type.II.Cells.1_LAM1','Monocytes_LAM1','Gamma.Delta.T.Cells_LAM1','Microglia_LAM1','Endothelial.Cells_LAM1','Endothelial.Cells.1_LAM1','Pulmonary.Alveolar.Type.II.Cells_LAM1','T.Cells_LAM1','Ependymal.Cells_LAM1','B.Cells_LAM1','Mast.Cells_LAM1','Endothelial.Cells.3_LAM1','NK.Cells_LAM1','Smooth.Muscle.Cells_LAM1','Erythroid.like.And.Erythroid.Precursor.Cells_LAM1'), ident.2 =c('Microglia_WT','Pulmonary.Alveolar.Type.II.Cells_WT','Endothelial.Cells.2_WT','Dendritic.Cells_WT','Endothelial.Cells.1_WT','NK.Cells_WT','B.Cells_WT','Endothelial.Cells.3_WT','Mast.Cells_WT','Dendritic.Cells.1_WT','T.Cells_WT','Monocytes_WT','Pulmonary.Alveolar.Type.II.Cells.1_WT','Endothelial.Cells_WT','Smooth.Muscle.Cells_WT','Erythroid.like.And.Erythroid.Precursor.Cells_WT','Gamma.Delta.T.Cells_WT','Ependymal.Cells_WT','Pulmonary.Alveolar.Type.II.Cells.2_WT'), verbose = FALSE)

indexsigbulk<-which(LAM_WT.markersbulk$p_val_adj<=0.05)
LAM_WT.markersbulk_sig<-LAM_WT.markersbulk[indexsigbulk,]

Up<-rownames(LAM_WT.markersbulk_sig[which(LAM_WT.markersbulk_sig$avg_log2FC>=0),])
Down<-rownames(LAM_WT.markersbulk_sig[which(LAM_WT.markersbulk_sig$avg_log2FC<0),])

#dbstouse <- c("PanglaoDB_Augmented_2021","CellMarker_Augmented_2021","Human_Gene_Atlas","Old_CMAP_up","Old_CMAP_down","LINCS_L1000_Chem_Pert_up","LINCS_L1000_Chem_Pert_down","LINCS_L1000_Ligand_Perturbations_down","LINCS_L1000_Ligand_Perturbations_up","KEGG_2021_Human","MSigDB_Hallmark_2020","GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015","ProteomicsDB_2020") #PanglaoDB_Augmented_2021
dbstouse <- c("Old_CMAP_up","Old_CMAP_down","LINCS_L1000_Chem_Pert_up","LINCS_L1000_Chem_Pert_down","LINCS_L1000_Ligand_Perturbations_down","LINCS_L1000_Ligand_Perturbations_up") #PanglaoDB_Augmented_2021
dbstousePaths<-c("KEGG_2021_Human","MSigDB_Hallmark_2020","GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015","ProteomicsDB_2020","PanglaoDB_Augmented_2021","CellMarker_Augmented_2021","Human_Gene_Atlas")

if(length(Up) > 200){
  Up<-Up[1:200]
}
if(length(Down) > 200){
  Down<-Down[1:200]
}

if (websiteLive) {
  enrichedUp <- enrichr(Up, dbstouse)
}

if (websiteLive) {
  enrichedDown <- enrichr(Down, dbstouse)
}
#names(enriched)
#celltype<-c(enriched[[1]]$Term[1],enriched[[2]]$Term[1])
#if (websiteLive) plotEnrich(enriched[["Old_CMAP_up"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CMAP ASTRO UP")+plotEnrich(enriched[["Old_CMAP_down"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CMAP ASTRO Down")
#if (websiteLive) plotEnrich(enriched[["KEGG_2021_Human"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "KEGG")+plotEnrich(enriched[["MSigDB_Hallmark_2020"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "MSigDB")+plotEnrich(enriched[["GO_Molecular_Function_2015"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "GO Molecular Function")

if (websiteLive) {
  enrichedPaths <- enrichr(c(Up,Down), dbstousePaths)
}
if (websiteLive) plotEnrich(enrichedPaths[["KEGG_2021_Human"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = paste("KEGG","LAM Bulk",sep=" "))#+plotEnrich(enriched[["MSigDB_Hallmark_2020"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "MSigDB")+plotEnrich(enriched[["GO_Molecular_Function_2015"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "GO Molecular Function")

pathsBULK<-enrichedPaths[["KEGG_2021_Human"]]

drugsUp<-enrichedUp[["Old_CMAP_down"]]$Term[which(enrichedUp[["Old_CMAP_down"]]$P.value <0.05)]
drugsUp<-gsub("-[0-9]+$","",drugsUp)

drugsDown<-enrichedDown[["Old_CMAP_up"]]$Term[which(enrichedDown[["Old_CMAP_up"]]$P.value <0.05)]
drugsDown<-gsub("-[0-9]+$","",drugsDown)
if (websiteLive) plotEnrich(enrichedUp[["Old_CMAP_down"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CMAP BULK UP")+plotEnrich(enrichedDown[["Old_CMAP_up"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CMAP BULK Down")

drugInfo<-read.delim("COVID/drug_repurposing_hub.txt")
AlldrugsCellIDs<-c(drugsUp,drugsDown)
foundds<-0
MOA<-c()
for(d in 1:length(AlldrugsCellIDs)){
  indexfoundd<-grep(AlldrugsCellIDs[d],drugInfo$pert_iname)
  if(!length(indexfoundd)==0){
    foundds<-foundds+1
    #MOA<-c(MOA,drugInfo$moa[indexfoundd])
    MOA<-rbind(MOA,cbind(drugInfo$pert_iname[indexfoundd],drugInfo$moa[indexfoundd],drugInfo$disease_area[indexfoundd]))
  }
}
indexesnull1<-which(MOA[,2] == "")
if(!length(indexesnull1)==0){
  MOA<-MOA[-indexesnull1,]
}
indexesnull2<-which(MOA[,3] == "")
if(!length(indexesnull2)==0){
  MOA<-MOA[-indexesnull2,]
}

# sort by Freq
tableMOAsBulk<-as.data.frame(table(MOA[,3]))
tableMOAsBulk <- tableMOAsBulk[order(tableMOAsBulk$Freq,decreasing = T),]
write.table(tableMOAsBulk,"tableMOAs_LAM_BULK.txt",quote = F,row.names = F,sep = "\t")

# if (websiteLive) {
#   enriched <- enrichr(Up[1:200], dbstouse)
# }
# 
# celltype<-c(enriched[[2]]$Term[1],enriched[[3]]$Term[1])
# if (websiteLive) plotEnrich(enriched[["Old_CMAP_up"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CMAP UP")+plotEnrich(enriched[["Old_CMAP_down"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CMAP Down")
# if (websiteLive) plotEnrich(enriched[["KEGG_2021_Human"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "KEGG")+plotEnrich(enriched[["MSigDB_Hallmark_2020"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "MSigDB")+plotEnrich(enriched[["GO_Molecular_Function_2015"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "GO Molecular Function")


#Concluding Results
plots <- VlnPlot(LAM1_WT.combined, features = c("FOS","HBB","SFTPC","APOE"), split.by = "orig.ident", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)

wrap_plots(plots = plots, ncol = 2,nrow = 2)

plots <- VlnPlot(LAM1_WT.combined, features = c("SCGB1A1","IGJ","SCGB3A1","CCL21"), split.by = "orig.ident", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)

wrap_plots(plots = plots, ncol = 2,nrow = 2)

##########################################################################################################


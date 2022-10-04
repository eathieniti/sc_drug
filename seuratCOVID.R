library(dplyr)
library(Seurat)
library(patchwork)
library(multtest)
library(metap)
library(ggplot2)
library(cowplot)
library(enrichR)
library(future)
setwd("/home/oulas/scRNA-Seq")

# Load the seurat object
choroid_plexus <- readRDS("COVID/COVID-19_brain_snRNA-seq_choroid_plexus_final_seurat_v3.2.3.rds")
parenchyma_cortex<- readRDS("COVID/COVID-19_brain_snRNA-seq_parenchyma_cortex_final_seurat_v3.2.3.rds")
unique(parenchyma_cortex$Sample)

#QC and selecting cells for further analysis
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#These have been run already!!!
#choroid_plexus[["percent.mt"]] <- PercentageFeatureSet(choroid_plexus, pattern = "^MT-")
#head(choroid_plexus@meta.data, 5)


# Visualize QC metrics as a violin plot
VlnPlot(choroid_plexus, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(choroid_plexus, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(choroid_plexus, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Filtering data
#choroid_plexus <- subset(choroid_plexus, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 75)


#Normalizing the data. Normalized values are stored in choroid_plexus[["RNA"]]@data. *
#These have been run already!!!
#choroid_plexus <- NormalizeData(choroid_plexus, normalization.method = "LogNormalize", scale.factor = 10000)


#Identification of highly variable features (feature selection) *
choroid_plexus <- FindVariableFeatures(choroid_plexus, selection.method = "vst", nfeatures = 2000)
parenchyma_cortex <- FindVariableFeatures(parenchyma_cortex, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(choroid_plexus), 10)


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(choroid_plexus)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#Scaling the data

all.genes <- rownames(choroid_plexus)
#These have been run already!!!
#choroid_plexus <- ScaleData(choroid_plexus, features = all.genes)


#Perform linear dimensional reduction
#These have been run already!!!
#choroid_plexus <- RunPCA(choroid_plexus, features = VariableFeatures(object = choroid_plexus))


# Examine and visualize PCA results a few different ways
print(choroid_plexus[["pca"]], dims = 1:5, nfeatures = 5)


VizDimLoadings(choroid_plexus, dims = 1:2, reduction = "pca")


DimPlot(choroid_plexus, reduction = "pca")


DimHeatmap(choroid_plexus, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(choroid_plexus, dims = 1:15, cells = 500, balanced = TRUE)


#Determine the dimensionality of the dataset

# NOTE: JackStraw process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
#choroid_plexus <- JackStraw(choroid_plexus, num.replicate = 100) #approx 3 mins
#choroid_plexus <- ScoreJackStraw(choroid_plexus, dims = 1:20)
#JackStrawPlot(choroid_plexus, dims = 1:15)
#These have been run already!!!
ElbowPlot(choroid_plexus)

#Cluster the cells
#These have been run already!!!
#choroid_plexus <- FindNeighbors(choroid_plexus, dims = 1:10)
#choroid_plexus <- FindClusters(choroid_plexus, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(choroid_plexus), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# umap-learn)
#These have been run already!!!
#choroid_plexus <- RunUMAP(choroid_plexus, dims = 1:10)
#choroid_plexus <- RunTSNE(choroid_plexus,dims = 1:10)

# note that you can set label = TRUE or use the LabelClusters function to help label
# individual clusters
DimPlot(choroid_plexus, reduction = "umap")
DimPlot(parenchyma_cortex, reduction = "umap")
DimPlot(choroid_plexus, reduction = "tsne")

DimPlot(choroid_plexus, reduction = "umap",group.by = "cellID",label = TRUE)
DimPlot(choroid_plexus, reduction = "umap", split.by = "Disease",label = TRUE)
DimPlot(choroid_plexus, reduction = "umap", split.by = "Biogroup",label = TRUE)

DimPlot(parenchyma_cortex, reduction = "umap", group.by = "cellID",label = TRUE)
DimPlot(parenchyma_cortex, reduction = "umap", split.by = "Disease",label = TRUE)
DimPlot(parenchyma_cortex, reduction = "umap", split.by = "Biogroup",label = TRUE)


#Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster Epithelial
clusterEp.markers <- FindMarkers(choroid_plexus, ident.1 = "Epithelial", min.pct = 0.25)
head(clusterEp.markers, n = 5)

unique(choroid_plexus$Biogroup)
# # find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster5.markers <- FindMarkers(choroid_plexus, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
# head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones - approx 4.5 mins with default test - approx 12 mins with MAST
#BiocManager::install("MAST")

choroid_plexus.markers <- FindAllMarkers(choroid_plexus, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)#,test.use = "MAST"
parenchyma_cortex.markers <- FindAllMarkers(parenchyma_cortex, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)#,test.use = "MAST"
parenchyma_cortex.markersAll <- FindAllMarkers(parenchyma_cortex, ident.1 = "Astrocyte", min.pct = 0.25)#,test.use = "MAST"
#saveRDS(choroid_plexus.markers, file = "choroid_plexus.markers.rds")
#choroid_plexus.markers<-readRDS(file = "choroid_plexus.markers.rds")

choroid_plexus.markerstop2<-choroid_plexus.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

choroid_plexus.markerstop100<-choroid_plexus.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC)



parenchyma_cortex.markerstop2<-parenchyma_cortex.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

parenchyma_cortex.markerstop100<-parenchyma_cortex.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC)


#Visualize markers choroid_plexus
VlnPlot(choroid_plexus, features = choroid_plexus.markerstop2$gene[choroid_plexus.markerstop2$cluster=="Epithelial"])


FeaturePlot(choroid_plexus, features = choroid_plexus.markerstop2$gene[choroid_plexus.markerstop2$cluster=="Epithelial"])


DoHeatmap(choroid_plexus, features = choroid_plexus.markerstop2$gene) + NoLegend()
DotPlot(choroid_plexus, features = choroid_plexus.markerstop2$gene, cols = c("blue", "red"), dot.scale = 8, split.by = "Biogroup") +
  RotatedAxis()



#Visualize markers parenchyma_cortex
VlnPlot(parenchyma_cortex, features = parenchyma_cortex.markerstop2$gene[parenchyma_cortex.markerstop2$cluster=="Astrocyte"])



FeaturePlot(parenchyma_cortex, features = parenchyma_cortex.markerstop2$gene[parenchyma_cortex.markerstop2$cluster=="Astrocyte"])


DoHeatmap(parenchyma_cortex, features = parenchyma_cortex.markerstop2$gene) + NoLegend()
DotPlot(parenchyma_cortex, features = parenchyma_cortex.markerstop2$gene, cols = c("blue", "red"), dot.scale = 8, split.by = "Biogroup") +
  RotatedAxis()


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
for(m in as.character(unique(choroid_plexus.markerstop100$cluster))){
  markers_to_search<-choroid_plexus.markerstop100$gene[choroid_plexus.markerstop100$cluster==m]
  if (websiteLive) {
    enriched <- enrichr(markers_to_search, dbstouse)
  }
  celltype<-c(enriched[["PanglaoDB_Augmented_2021"]]$Term[1],enriched[["CellMarker_Augmented_2021"]]$Term[1])
  print(paste("Cluster",m,celltype,sep = " "))
  names(enriched)<-paste(names(enriched),m,sep = "_")
  listofresults<-c(listofresults,enriched)
  annotatedclusters<-c(annotatedclusters,toString(celltype[2]))
}
if (websiteLive) plotEnrich(listofresults[["PanglaoDB_Augmented_2021_Epithelial"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "PanglaoDB")+plotEnrich(listofresults[["CellMarker_Augmented_2021_Epithelial"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CellMarker")
if (websiteLive) plotEnrich(listofresults[["PanglaoDB_Augmented_2021_Mesenchymal"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "PanglaoDB")+plotEnrich(listofresults[["CellMarker_Augmented_2021_Mesenchymal"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CellMarker")
if (websiteLive) plotEnrich(listofresults[["PanglaoDB_Augmented_2021_Neural"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "PanglaoDB")+plotEnrich(listofresults[["CellMarker_Augmented_2021_Neural"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CellMarker")
if (websiteLive) plotEnrich(listofresults[["PanglaoDB_Augmented_2021_Macrophage"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "PanglaoDB")+plotEnrich(listofresults[["CellMarker_Augmented_2021_Macrophage"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CellMarker")
if (websiteLive) plotEnrich(listofresults[["PanglaoDB_Augmented_2021_Ependymal"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "PanglaoDB")+plotEnrich(listofresults[["CellMarker_Augmented_2021_Ependymal"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CellMarker")
if (websiteLive) plotEnrich(listofresults[["PanglaoDB_Augmented_2021_Glial"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "PanglaoDB")+plotEnrich(listofresults[["CellMarker_Augmented_2021_Glial"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CellMarker")
if (websiteLive) plotEnrich(listofresults[["PanglaoDB_Augmented_2021_Endothelial"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "PanglaoDB")+plotEnrich(listofresults[["CellMarker_Augmented_2021_Endothelial"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CellMarker")

print(annotatedclusters)


parenchyma_cortex$Biogroup.cellID <- paste(parenchyma_cortex$Biogroup, parenchyma_cortex$cellID, sep = "_")
choroid_plexus$Biogroup.cellID <- paste(choroid_plexus$Biogroup, choroid_plexus$cellID, sep = "_")
#choroid_plexus_parenchyma_cortex.combined$celltype <- Idents(choroid_plexus_parenchyma_cortex.combined)
Idents(parenchyma_cortex) <- "Biogroup.cellID"
Idents(choroid_plexus) <- "Biogroup.cellID"
unique(Idents(parenchyma_cortex))
unique(Idents(choroid_plexus))

for(cellIDs in as.character(unique(parenchyma_cortex$cellID))){
  
}

#Already run I think???
parenchyma_cortex <- NormalizeData(parenchyma_cortex, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(parenchyma_cortex)
parenchyma_cortex <- ScaleData(parenchyma_cortex, features = all.genes)
DefaultAssay(parenchyma_cortex) <- "integrated"
#Already run I think???
choroid_plexus <- NormalizeData(choroid_plexus, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(choroid_plexus)
choroid_plexus <- ScaleData(choroid_plexus, features = all.genes)
DefaultAssay(choroid_plexus) <- "integrated"


cellIDs<-"Astrocyte"
cellIDs<-"Endothelial"
Covid<-parenchyma_cortex$Biogroup.cellID[grep(paste("Case_",cellIDs,sep=""),parenchyma_cortex$Biogroup.cellID)]
Control<-parenchyma_cortex$Biogroup.cellID[grep(paste("Control_",cellIDs,sep=""),parenchyma_cortex$Biogroup.cellID)]

Covid<-choroid_plexus$Biogroup.cellID[grep(paste("Case_",cellIDs,sep=""),choroid_plexus$Biogroup.cellID)]
Control<-choroid_plexus$Biogroup.cellID[grep(paste("Control_",cellIDs,sep=""),choroid_plexus$Biogroup.cellID)]

DefaultAssay(parenchyma_cortex) <- "RNA"
parenchyma_cortex.markersCellIDs <- FindMarkers(parenchyma_cortex, ident.1 = Control, ident.2 = Covid, verbose = FALSE)
choroid_plexus.markersCellIDs <- FindMarkers(choroid_plexus, ident.1 = Control, ident.2 = Covid, verbose = FALSE)

astros<-Idents(parenchyma_cortex)[grep("Astrocyte",parenchyma_cortex$cellID)]
parenchyma_cortexAstro<-subset(parenchyma_cortex, idents = astros)
DefaultAssay(parenchyma_cortexAstro) <- "RNA"#is this needed???

endothelial<-Idents(choroid_plexus)[grep("Endothelial",choroid_plexus$cellID)]
choroid_plexusEndo<-subset(choroid_plexus, idents = endothelial)
DefaultAssay(choroid_plexusEndo) <- "RNA"#is this needed???

indexsigclCellIDs<-which(parenchyma_cortex.markersCellIDs$p_val_adj<=0.05)
parenchyma_cortex.markersCellIDs_sig<-parenchyma_cortex.markersCellIDs[indexsigclCellIDs,]
#which(rownames(parenchyma_cortex.markersCellIDs_sig)=="FTH1")
Up<-rownames(parenchyma_cortex.markersCellIDs_sig[which(parenchyma_cortex.markersCellIDs_sig$avg_log2FC>=0),])
Down<-rownames(parenchyma_cortex.markersCellIDs_sig[which(parenchyma_cortex.markersCellIDs_sig$avg_log2FC<0),])

indexsigclCellIDs<-which(choroid_plexus.markersCellIDs$p_val_adj<=0.05)
choroid_plexus.markersCellIDs_sig<-choroid_plexus.markersCellIDs[indexsigclCellIDs,]
#which(rownames(parenchyma_cortex.markersCellIDs_sig)=="FTH1")
Up<-rownames(choroid_plexus.markersCellIDs_sig[which(choroid_plexus.markersCellIDs_sig$avg_log2FC>=0),])
Down<-rownames(choroid_plexus.markersCellIDs_sig[which(choroid_plexus.markersCellIDs_sig$avg_log2FC<0),])

choroid_plexusEndo <- ScaleData(choroid_plexusEndo, features = all.genes)
DoHeatmap(choroid_plexusEndo, features = c(Up,Down),raster=F,size = 4)+
  theme(axis.text.y = element_text(size = 8))

length(Up)+length(Down)
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

if (websiteLive) plotEnrich(enrichedPaths[["KEGG_2021_Human"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = paste("KEGG",cellIDs,sep=" "))
paths<-enrichedPaths[["KEGG_2021_Human"]]

drugsUp<-enrichedUp[["Old_CMAP_down"]]$Term[which(enrichedUp[["Old_CMAP_down"]]$P.value <0.05)]
drugsUp<-gsub("-[0-9]+$","",drugsUp)

drugsDown<-enrichedDown[["Old_CMAP_up"]]$Term[which(enrichedDown[["Old_CMAP_up"]]$P.value <0.05)]
drugsDown<-gsub("-[0-9]+$","",drugsDown)

if (websiteLive) plotEnrich(enrichedUp[["Old_CMAP_down"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = paste("CMAP", cellIDs, "UP",sep=" "))+plotEnrich(enrichedDown[["Old_CMAP_up"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = paste("CMAP", cellIDs, "Down",sep=" "))

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
write.table(tableMOAs,paste("tableMOAs_",cellIDs,".txt",sep=""),quote = F,row.names = F,sep = "\t")


#Bulk
CVBulk<-parenchyma_cortex$Biogroup.cellID[grep("Case",parenchyma_cortex$Biogroup.cellID)]
CTBulk<-parenchyma_cortex$Biogroup.cellID[grep("Control",parenchyma_cortex$Biogroup.cellID)]
plan()
plan("multiprocess", workers = 4)

#need to add this to run 
options(future.globals.maxSize = 8000 * 1024^2)
parenchyma_cortex.markersBulk <- FindMarkers(parenchyma_cortex, ident.1 = CTBulk, ident.2 = CVBulk, verbose = FALSE)
indexsigclBulk<-which(parenchyma_cortex.markersBulk$p_val_adj<=0.05)
parenchyma_cortex.markersBulk_sig<-parenchyma_cortex.markersBulk[indexsigclBulk,]

Up<-rownames(parenchyma_cortex.markersBulk_sig[which(parenchyma_cortex.markersBulk_sig$avg_log2FC>=0),])
Down<-rownames(parenchyma_cortex.markersBulk_sig[which(parenchyma_cortex.markersBulk_sig$avg_log2FC<0),])
length(Up)+length(Down)
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

if (websiteLive) {
  enrichedPaths <- enrichr(c(Up,Down), dbstousePaths)
}

#names(enriched)
#celltype<-c(enriched[[1]]$Term[1],enriched[[2]]$Term[1])
#if (websiteLive) plotEnrich(enriched[["Old_CMAP_up"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CMAP BULK UP")+plotEnrich(enriched[["Old_CMAP_down"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CMAP BULK Down")
#if (websiteLive) plotEnrich(enriched[["KEGG_2021_Human"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "KEGG")+plotEnrich(enriched[["MSigDB_Hallmark_2020"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "MSigDB")+plotEnrich(enriched[["GO_Molecular_Function_2015"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "GO Molecular Function")

paths<-enrichedPaths[["KEGG_2021_Human"]]$Term
if (websiteLive) plotEnrich(enrichedPaths[["KEGG_2021_Human"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = paste("KEGG","Bulk",sep=" "))#+plotEnrich(enriched[["MSigDB_Hallmark_2020"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "MSigDB")+plotEnrich(enriched[["GO_Molecular_Function_2015"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "GO Molecular Function")

drugsUp<-enrichedUp[["Old_CMAP_down"]]$Term[which(enrichedUp[["Old_CMAP_down"]]$P.value <0.05)]
drugsUp<-gsub("-[0-9]+$","",drugsUp)

drugsDown<-enrichedDown[["Old_CMAP_up"]]$Term[which(enrichedDown[["Old_CMAP_up"]]$P.value <0.05)]
drugsDown<-gsub("-[0-9]+$","",drugsDown)

if (websiteLive) plotEnrich(enrichedUp[["Old_CMAP_down"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CMAP BULK UP")+plotEnrich(enrichedDown[["Old_CMAP_up"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CMAP BULK Down")

drugInfo<-read.delim("COVID/drug_repurposing_hub.txt")
AlldrugsBulk<-c(drugsUp,drugsDown)
foundds<-0
MOABulk<-c()
for(d in 1:length(AlldrugsBulk)){
  indexfoundd<-which(startsWith(AlldrugsBulk[d],drugInfo$pert_iname)==TRUE)
  
  if(!length(indexfoundd)==0){
    foundds<-foundds+1
    #MOABulk<-c(MOABulk,drugInfo$moa[indexfoundd])
    #if(!drugInfo$moa[indexfoundd]==""){# !drugInfo$disease_area[indexfoundd] == ""){
      MOABulk<-rbind(MOABulk,cbind(drugInfo$pert_iname[indexfoundd],drugInfo$moa[indexfoundd],drugInfo$disease_area[indexfoundd]))
    #}
  }
}
indexesnull1<-which(MOABulk[,2] == "")
if(!length(indexesnull1)==0){
  MOABulk<-MOABulk[-indexesnull1,]
}

indexesnull2<-which(MOABulk[,3] == "")
if(!length(indexesnull2)==0){
  MOABulk<-MOABulk[-indexesnull2,]
}


# sort by Freq
tableMOAsBulk<-as.data.frame(table(MOABulk[,3]))
tableMOAsBulk <- tableMOAsBulk[order(tableMOAsBulk$Freq,decreasing = T),]
write.table(tableMOAsBulk,"tableMOAs_BULK.txt",quote = F,row.names = F,sep = "\t")

#Change labels
new.cluster.ids <- make.names(annotatedclusters,unique = T)
names(new.cluster.ids) <- levels(choroid_plexus)
choroid_plexus <- RenameIdents(choroid_plexus, new.cluster.ids)
unique(Idents(choroid_plexus))
DimPlot(choroid_plexus, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
############################################################################################################################



#Run Integrated analysis#############################################################################
choroid_plexus <- readRDS("COVID/COVID-19_brain_snRNA-seq_choroid_plexus_final_seurat_v3.2.3.rds")
parenchyma_cortex<- readRDS("COVID/COVID-19_brain_snRNA-seq_parenchyma_cortex_final_seurat_v3.2.3.rds")

#choroid_plexus[["percent.mt"]] <- PercentageFeatureSet(choroid_plexus, pattern = "^MT-")
#parenchyma_cortex[["percent.mt"]] <- PercentageFeatureSet(parenchyma_cortex, pattern = "^MT-")
choroid_plexus_parenchyma_cortex.list<-c(choroid_plexus,parenchyma_cortex)


# choroid_plexus_parenchyma_cortex.list <- lapply(X = choroid_plexus_parenchyma_cortex.list, FUN = function(x) {
#   #Filtering data
#   #x <- subset(x, subset = nFeature_RNA > 100 & nFeature_RNA < 7000 & percent.mt < 85)
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# })


features <- SelectIntegrationFeatures(object.list = choroid_plexus_parenchyma_cortex.list)
choroid_plexus_parenchyma_cortex.anchors <- FindIntegrationAnchors(object.list = choroid_plexus_parenchyma_cortex.list, dims = 1:20)#anchor.features = features
choroid_plexus_parenchyma_cortex.combined <- IntegrateData(anchorset = choroid_plexus_parenchyma_cortex.anchors, dims = 1:20)
unique(choroid_plexus_parenchyma_cortex.combined$Sample)
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(choroid_plexus_parenchyma_cortex.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
choroid_plexus_parenchyma_cortex.combined <- ScaleData(choroid_plexus_parenchyma_cortex.combined, verbose = FALSE)
choroid_plexus_parenchyma_cortex.combined <- RunPCA(choroid_plexus_parenchyma_cortex.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
choroid_plexus_parenchyma_cortex.combined <- RunUMAP(choroid_plexus_parenchyma_cortex.combined, reduction = "pca", dims = 1:20)
choroid_plexus_parenchyma_cortex.combined <- RunTSNE(choroid_plexus_parenchyma_cortex.combined,reduction = "pca",dims = 1:10)
#Louvein method or improved Leiden  ??
choroid_plexus_parenchyma_cortex.combined <- FindNeighbors(choroid_plexus_parenchyma_cortex.combined, reduction = "pca", dims = 1:20)
choroid_plexus_parenchyma_cortex.combined <- FindClusters(choroid_plexus_parenchyma_cortex.combined, resolution  = 0.3)
#Can view object
choroid_plexus_parenchyma_cortex.combined[[]]
#choroid_plexus_parenchyma_cortex.combined[['groups']] <- sample(x = c('g1', 'g2'), size = ncol(x = choroid_plexus_parenchyma_cortex.combined), replace = TRUE)
p1combined <- DimPlot(choroid_plexus_parenchyma_cortex.combined, reduction = "umap", group.by = "orig.ident")
p2combined <- DimPlot(choroid_plexus_parenchyma_cortex.combined, reduction = "umap", label = TRUE)
p1combined+p2combined



DimPlot(choroid_plexus_parenchyma_cortex.combined, reduction = "umap", split.by = "Disease")
DimPlot(choroid_plexus_parenchyma_cortex.combined, reduction = "tsne", split.by = "Disease")


#Identify conserved cell type markers
# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(choroid_plexus_parenchyma_cortex.combined) <- "RNA"



choroid_plexus_parenchyma_cortex.markers <- FindAllMarkers(choroid_plexus_parenchyma_cortex.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
choroid_plexus_parenchyma_cortex.markerstop1<-choroid_plexus_parenchyma_cortex.markers %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC) #100 gives best results

choroid_plexus_parenchyma_cortex.markerstop2<-choroid_plexus_parenchyma_cortex.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC) #100 gives best results

choroid_plexus_parenchyma_cortex.markerstop100<-choroid_plexus_parenchyma_cortex.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC) #100  gives best results

choroid_plexus_parenchyma_cortex.markerstop200<-choroid_plexus_parenchyma_cortex.markers %>%
  group_by(cluster) %>%
  slice_max(n = 200, order_by = avg_log2FC) #200 gives second best results


FeaturePlot(choroid_plexus_parenchyma_cortex.combined, features = choroid_plexus_parenchyma_cortex.markerstop1$gene[c(1,2,3,4,5,6,9,10,21)])


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
for(m in as.numeric(as.character(unique(choroid_plexus_parenchyma_cortex.markerstop100$cluster)))){
  markers_to_search<-choroid_plexus_parenchyma_cortex.markerstop100$gene[choroid_plexus_parenchyma_cortex.markerstop100$cluster==m]
  if (websiteLive) {
    enriched <- enrichr(markers_to_search, dbstouse)
  }
  celltype<-c(enriched[[1]]$Term[1],enriched[[2]]$Term[1])
  print(paste("Cluster",m,celltype,sep = " "))
  annotatedclusters<-c(annotatedclusters,toString(celltype[2]))
  names(enriched)<-paste(names(enriched),m,sep = "_")
  listofresults<-c(listofresults,enriched)
}

if (websiteLive) plotEnrich(listofresults[["PanglaoDB_Augmented_2021_9"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "PanglaoDB")+plotEnrich(listofresults[["CellMarker_Augmented_2021_9"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CellMarker")

new.cluster.ids <- make.names(annotatedclusters,unique = T)
names(new.cluster.ids) <- levels(choroid_plexus_parenchyma_cortex.combined)
choroid_plexus_parenchyma_cortex.combined <- RenameIdents(choroid_plexus_parenchyma_cortex.combined, new.cluster.ids)


DimPlot(choroid_plexus_parenchyma_cortex.combined, reduction = "umap",label = TRUE)

markers.to.plot <- unique(choroid_plexus_parenchyma_cortex.markerstop2$gene) 
DotPlot(choroid_plexus_parenchyma_cortex.combined, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "orig.ident") +
  RotatedAxis()


#Change the Labels
choroid_plexus_parenchyma_cortex.combined$celltype.orig.ident <- paste(Idents(choroid_plexus_parenchyma_cortex.combined), choroid_plexus_parenchyma_cortex.combined$orig.ident, sep = "_")
choroid_plexus_parenchyma_cortex.combined$celltype <- Idents(choroid_plexus_parenchyma_cortex.combined)
Idents(choroid_plexus_parenchyma_cortex.combined) <- "celltype.orig.ident"
unique(Idents(choroid_plexus_parenchyma_cortex.combined))

CVEp<-Idents(choroid_plexus_parenchyma_cortex.combined)[grep("Epithelial.+_cv_",Idents(choroid_plexus_parenchyma_cortex.combined))]
CTEp<-Idents(choroid_plexus_parenchyma_cortex.combined)[grep("Epithelial.+_ct_",Idents(choroid_plexus_parenchyma_cortex.combined))]
choroid_plexus_parenchyma_cortex.markersSM <- FindMarkers(choroid_plexus_parenchyma_cortex.combined, ident.1 = CVEp, ident.2 = CTEp, verbose = FALSE)
indexsigclSM<-which(choroid_plexus_parenchyma_cortex.markersSM$p_val_adj<=0.05)
choroid_plexus_parenchyma_cortex.markersSM_sig<-choroid_plexus_parenchyma_cortex.markersSM[indexsigclSM,]

Up<-rownames(choroid_plexus_parenchyma_cortex.markersSM_sig[which(choroid_plexus_parenchyma_cortex.markersSM_sig$avg_log2FC>=0),])
Down<-rownames(choroid_plexus_parenchyma_cortex.markersSM_sig[which(choroid_plexus_parenchyma_cortex.markersSM_sig$avg_log2FC<0),])

dbstouse <- c("PanglaoDB_Augmented_2021","CellMarker_Augmented_2021","Human_Gene_Atlas","Old_CMAP_up","Old_CMAP_down","LINCS_L1000_Chem_Pert_up","LINCS_L1000_Chem_Pert_down","LINCS_L1000_Ligand_Perturbations_down","LINCS_L1000_Ligand_Perturbations_up","KEGG_2021_Human","MSigDB_Hallmark_2020","GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015","ProteomicsDB_2020") #PanglaoDB_Augmented_2021

if (websiteLive) {
  enriched <- enrichr(Up[1:200], dbstouse)
}
names(enriched)
celltype<-c(enriched[[1]]$Term[1],enriched[[2]]$Term[1])
if (websiteLive) plotEnrich(enriched[["Old_CMAP_up"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CMAP UP")+plotEnrich(enriched[["Old_CMAP_down"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "CMAP Down")
if (websiteLive) plotEnrich(enriched[["KEGG_2021_Human"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "KEGG")+plotEnrich(enriched[["MSigDB_Hallmark_2020"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "MSigDB")+plotEnrich(enriched[["GO_Molecular_Function_2015"]], showTerms = 40, numChar = 40, y = "Count", orderBy = "P.value",title = "GO Molecular Function")

drugsUp<-enriched[["Old_CMAP_up"]]$Term[which(enriched[["Old_CMAP_up"]]$Adjusted.P.value <0.05)]
drugsUp<-gsub("-[0-9]+$","",drugsUp)

drugInfo<-read.delim("COVID/drug_repurposing_hub.txt")

foundds<-0
MOA<-c()
for(d in 1:length(drugsUp)){
  indexfoundd<-grep(drugsUp[d],drugInfo$pert_iname)
  if(!length(indexfoundd)==0){
    foundds<-foundds+1
    MOA<-c(MOA,drugInfo$moa[indexfoundd])
  }
}


# sort by Freq
tableMOAs<-as.data.frame(table(MOA))
tableMOAs <- tableMOAs[order(tableMOAs$Freq,decreasing = T),]

#Subsetting cluster
#Identify differential expressed genes across conditions
choroid_plexus_parenchyma_cortex.combinedclSM<-subset(choroid_plexus_parenchyma_cortex.combined, idents = c("Smooth.Muscle.Cells_choroid_plexus","Smooth.Muscle.Cells_parenchyma_cortex"))
#Idents(choroid_plexus_parenchyma_cortex.combinedcl9)<-"choroid_plexus"
choroid_plexus_parenchyma_cortex.combinedclSM<-ScaleData(choroid_plexus_parenchyma_cortex.combinedclSM, features = rownames(choroid_plexus_parenchyma_cortex.markersSM_sig))

theme_set(theme_cowplot())
avg.clSM <- as.data.frame(log1p(AverageExpression(choroid_plexus_parenchyma_cortex.combinedclSM, verbose = FALSE)$RNA))
avg.clSM$gene <- rownames(avg.clSM)
genes.to.label = rownames(choroid_plexus_parenchyma_cortex.markersSM_sig)[1:10]
avg.clSM$gene[which(!avg.clSM$gene%in%genes.to.label)]<-""

p1 <- ggplot(avg.clSM, aes(`Smooth.Muscle.Cells_choroid_plexus`, `Smooth.Muscle.Cells_parenchyma_cortex`,label=gene)) + geom_point() + ggtitle("choroid_plexus signature genes")+geom_text()
p1

DoHeatmap(choroid_plexus_parenchyma_cortex.combinedclSM, features = c(Up[1:50],Down),raster=F,size = 4)+
  theme(axis.text.y = element_text(size = 8))

#Simulation Bulk RNA-seq
CVCells<-unique(Idents(choroid_plexus_parenchyma_cortex.combined))[grep("_cv_",unique(Idents(choroid_plexus_parenchyma_cortex.combined)))]
CTCells<-unique(Idents(choroid_plexus_parenchyma_cortex.combined))[grep("_ct_",unique(Idents(choroid_plexus_parenchyma_cortex.combined)))]
choroid_plexus_parenchyma_cortex.markersbulk <- FindMarkers(choroid_plexus_parenchyma_cortex.combined, ident.1 = CVCells,ident.2 = CTCells, verbose = FALSE)
#                                    c('Dendritic.Cells_choroid_plexus','Dendritic.Cells.1_choroid_plexus','Pulmonary.Alveolar.Type.II.Cells.2_choroid_plexus','Pulmonary.Alveolar.Type.II.Cells.1_choroid_plexus','Monocytes_choroid_plexus','Gamma.Delta.T.Cells_choroid_plexus','Microglia_choroid_plexus','Endothelial.Cells_choroid_plexus','Endothelial.Cells.1_choroid_plexus','Pulmonary.Alveolar.Type.II.Cells_choroid_plexus','T.Cells_choroid_plexus','Ependymal.Cells_choroid_plexus','B.Cells_choroid_plexus','Mast.Cells_choroid_plexus','Endothelial.Cells.3_choroid_plexus','NK.Cells_choroid_plexus','Smooth.Muscle.Cells_choroid_plexus','Erythroid.like.And.Erythroid.Precursor.Cells_choroid_plexus'), ident.2 =c('Microglia_parenchyma_cortex','Pulmonary.Alveolar.Type.II.Cells_parenchyma_cortex','Endothelial.Cells.2_parenchyma_cortex','Dendritic.Cells_parenchyma_cortex','Endothelial.Cells.1_parenchyma_cortex','NK.Cells_parenchyma_cortex','B.Cells_parenchyma_cortex','Endothelial.Cells.3_parenchyma_cortex','Mast.Cells_parenchyma_cortex','Dendritic.Cells.1_parenchyma_cortex','T.Cells_parenchyma_cortex','Monocytes_parenchyma_cortex','Pulmonary.Alveolar.Type.II.Cells.1_parenchyma_cortex','Endothelial.Cells_parenchyma_cortex','Smooth.Muscle.Cells_parenchyma_cortex','Erythroid.like.And.Erythroid.Precursor.Cells_parenchyma_cortex','Gamma.Delta.T.Cells_parenchyma_cortex','Ependymal.Cells_parenchyma_cortex','Pulmonary.Alveolar.Type.II.Cells.2_parenchyma_cortex'), verbose = FALSE)

indexsigbulk<-which(choroid_plexus_parenchyma_cortex.markersbulk$p_val_adj<=0.05)
choroid_plexus_parenchyma_cortex.markersbulk_sig<-choroid_plexus_parenchyma_cortex.markersbulk[indexsigbulk,]


Up<-rownames(choroid_plexus_parenchyma_cortex.markersbulk_sig[which(choroid_plexus_parenchyma_cortex.markersbulk_sig$avg_log2FC>=0),])

dbstouse <- c("PanglaoDB_Augmented_2021","CellMarker_Augmented_2021","Human_Gene_Atlas","Old_CMAP_up","Old_CMAP_down","LINCS_L1000_Chem_Pert_up","LINCS_L1000_Chem_Pert_down","LINCS_L1000_Ligand_Perturbations_down","LINCS_L1000_Ligand_Perturbations_up","KEGG_2021_Human","MSigDB_Hallmark_2020","GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015","ProteomicsDB_2020") #PanglaoDB_Augmented_2021

if (websiteLive) {
  enriched <- enrichr(Up[1:200], dbstouse)
}

celltype<-c(enriched[[2]]$Term[1],enr
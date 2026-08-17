#@@@@@@@@@@@@@@@@@@@@@@@
# ==============================================================================
# SPATIAL TRANSCRIPTOMICS MULTI-SAMPLE ANALYSIS PIPELINE
# Adjusted & Corrected Production Script
# ==============================================================================

rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringi)
library(SpotClean)
library(S4Vectors)

s_dir <- "/home/mekalav/Spatial/UABBCM/Input/spaceranger_outputs/"
out_dir <- "/home/mekalav/Spatial/output/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# STEP 1: Load Raw 10X Spatial Samples
# ------------------------------------------------------------------------------
b01 <- Load10X_Spatial(paste0(s_dir, "B01/outs"), slice = "B01"); b01$sample_id <- "B01"
b02 <- Load10X_Spatial(paste0(s_dir, "B02/outs"), slice = "B02"); b02$sample_id <- "B02"
b03 <- Load10X_Spatial(paste0(s_dir, "B03/outs"), slice = "B03"); b03$sample_id <- "B03"
b04 <- Load10X_Spatial(paste0(s_dir, "B04/outs"), slice = "B04"); b04$sample_id <- "B04"
b05 <- Load10X_Spatial(paste0(s_dir, "B05/outs"), slice = "B05"); b05$sample_id <- "B05"
b06 <- Load10X_Spatial(paste0(s_dir, "B06/outs"), slice = "B06"); b06$sample_id <- "B06"
b07 <- Load10X_Spatial(paste0(s_dir, "B07/outs"), slice = "B07"); b07$sample_id <- "B07"
b08 <- Load10X_Spatial(paste0(s_dir, "B08/outs"), slice = "B08"); b08$sample_id <- "B08"
b09 <- Load10X_Spatial(paste0(s_dir, "B09/outs"), slice = "B09"); b09$sample_id <- "B09"
b10 <- Load10X_Spatial(paste0(s_dir, "B10/outs"), slice = "B10"); b10$sample_id <- "B10"
b11 <- Load10X_Spatial(paste0(s_dir, "B11/outs"), slice = "B11"); b11$sample_id <- "B11"

w01 <- Load10X_Spatial(paste0(s_dir, "W01/outs"), slice = "W01"); w01$sample_id <- "W01"
w02 <- Load10X_Spatial(paste0(s_dir, "W02/outs"), slice = "W02"); w02$sample_id <- "W02"
w03 <- Load10X_Spatial(paste0(s_dir, "W03/outs"), slice = "W03"); w03$sample_id <- "W03"
w04 <- Load10X_Spatial(paste0(s_dir, "W04/outs"), slice = "W04"); w04$sample_id <- "W04"
w05 <- Load10X_Spatial(paste0(s_dir, "W05/outs"), slice = "W05"); w05$sample_id <- "W05"
w06 <- Load10X_Spatial(paste0(s_dir, "W06/outs"), slice = "W06"); w06$sample_id <- "W06"
w07 <- Load10X_Spatial(paste0(s_dir, "W07/outs"), slice = "W07"); w07$sample_id <- "W07"
w08 <- Load10X_Spatial(paste0(s_dir, "W08/outs"), slice = "W08"); w08$sample_id <- "W08"
w09 <- Load10X_Spatial(paste0(s_dir, "W09/outs"), slice = "W09"); w09$sample_id <- "W09"
w10 <- Load10X_Spatial(paste0(s_dir, "W10/outs"), slice = "W10"); w10$sample_id <- "W10"
w11 <- Load10X_Spatial(paste0(s_dir, "W11/outs"), slice = "W11"); w11$sample_id <- "W11"

# Combine into a single macro-object
combined_seurat <- merge(
  x = b01, 
  y = c(b02, b03, b04, b05, b06, b07, b08, b09, b10, b11, w01, w02, w03, w04, w05, w06, w07, w08, w09, w10, w11),
  add.cell.ids = c("B01","B02","B03","B04","B05","B06","B07","B08","B09","B10","B11","W01","W02","W03","W04","W05","W06","W07","W08","W09","W10","W11"),
  project = "TNBC_Project")
print(combined_seurat)

rm(b01,b02,b03,b04,b05,b06,b07,b08,b09,b10,b11,w01,w02,w03,w04,w05,w06,w07,w08,w09,w10,w11); gc()

# ------------------------------------------------------------------------------
# STEP 2: Demographic Metadata Mapping & QC Filtering
# ------------------------------------------------------------------------------
metadata <- read.csv("/home/mekalav/Spatial/UABBCM/Input/ADImetadata.csv", stringsAsFactors = FALSE)
head(metadata)
columns_to_add <- setdiff(colnames(metadata), "De_ID")
for (col in columns_to_add) {
  combined_seurat[[col]] <- metadata[[col]][match(combined_seurat$sample_id, metadata$De_ID)]}

head(combined_seurat@meta.data)

# Compute mitochondrial percentage
combined_seurat[["percent.mt"]] <- PercentageFeatureSet(combined_seurat, pattern = "^MT-")

cat("Total spots before filtering:", ncol(combined_seurat), "\n")
combined_seurat <- subset(
  combined_seurat,
  subset = nFeature_Spatial > 300 & nFeature_Spatial < 8000 & percent.mt < 20)
cat("Total spots after filtering:", ncol(combined_seurat), "\n")

# ------------------------------------------------------------------------------
# STEP 3: Generate Base QC Metrics Plots
# ------------------------------------------------------------------------------
cohort_summary <- combined_seurat@meta.data %>%
  group_by(sample_id) %>%
  summarise(
    Total_Spots  = n(),
    Median_Genes = median(nFeature_Spatial),
    Median_UMIs  = median(nCount_Spatial),
    Race         = unique(Race)
  )

spot_plot <- ggplot(cohort_summary, aes(x = sample_id, y = Total_Spots, fill = Race)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  scale_fill_manual(values = c("Black" = "#E64B35FF", "White" = "#4DBBD5FF")) +
  theme_classic(base_size = 12) +
  labs(title = "Cohort Overview: Spot-Level Spatial Throughput", x = "Visium Sample ID", y = "Total Captured Spots") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(filename = paste0(out_dir, "P1_Cohort_Spot_Level_Total_Barplot.pdf"), plot = spot_plot, width = 10, height = 5)

p2 <- VlnPlot(object = combined_seurat, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), group.by = "sample_id", pt.size = 0, ncol = 1)
ggsave(filename = paste0(out_dir, "P2_vlnplot.pdf"), plot = p2, width = 8, height = 10)
# ------------------------------------------------------------------------------
# STEP 4: High-Speed Layer Normalization (SCTransform)
# ------------------------------------------------------------------------------
# Split raw spatial layers by sample to normalize each batch independently
combined_seurat[["Spatial"]] <- split(combined_seurat[["Spatial"]], f = combined_seurat$sample_id)
View(combined_seurat@meta.data)
head(combined_seurat@meta.data)

combined_seurat <- SCTransform(
  object = combined_seurat, 
  assay = "Spatial", 
  vars.to.regress = "percent.mt", 
  variable.features.n = 3000, 
  verbose = TRUE)

saveRDS(combined_seurat, paste0(out_dir, "integrated_data.rds"))




# ------------------------------------------------------------------------------
# STEP 5: Base Dimensionality Reduction & Layer Integration (RPCA)
# ------------------------------------------------------------------------------
# Join layers temporarily so PCA calculates variations across all 51,647 spots at once
combined_seurat[["SCT"]] <- JoinLayers(combined_seurat[["SCT"]])
combined_seurat <- RunPCA(combined_seurat, assay = "SCT", verbose = FALSE)

# Re-split the normalized assay layers so the integration engine knows batch borders
combined_seurat[["SCT"]] <- split(combined_seurat[["SCT"]], f = combined_seurat$sample_id)

# Run modern Seurat v5 integration
combined_seurat <- IntegrateLayers(
  object = combined_seurat, 
  method = RPCAIntegration, 
  orig.reduction = "pca",
  new.reduction = "integrated.rpca", 
  verbose = TRUE)

Images(combined_seurat)
head(GetTissueCoordinates(combined_seurat))

# ------------------------------------------------------------------------------
# STEP 6: Community Graph Clustering & UMAP Projections
# ------------------------------------------------------------------------------
combined_seurat <- FindNeighbors(combined_seurat, reduction = "integrated.rpca", dims = 1:30)

library(clustree)
combined_seurat <- FindClusters(
  object = combined_seurat, 
  resolution = seq(0, 1.5, by = 0.1))

pdf("/home/mekalav/Spatial/output/Clustree_Resolution_Check.pdf", width = 10, height = 12)
clustree(combined_seurat, prefix = "SCT_snn_res.")
dev.off()

combined_seurat <- FindClusters(combined_seurat, resolution = 0.5) 
combined_seurat <- RunUMAP(combined_seurat, reduction = "integrated.rpca", dims = 1:30)

saveRDS(combined_seurat, paste0(out_dir, "integrated_data1.rds"))

# ------------------------------------------------------------------------------
# STEP 7: Export Multi-Panel Evaluation Visualizations
# ------------------------------------------------------------------------------
# Reset identities and sync generic name tags
Idents(combined_seurat) <- "seurat_clusters"
combined_seurat$orig.ident <- combined_seurat$sample_id

p1_umap <- DimPlot(combined_seurat, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 4, pt.size = 0.1) + 
  theme_classic(base_size = 11) + labs(title = "A. Spatial Niche Clusters") + NoLegend()

p2_umap <- DimPlot(combined_seurat, reduction = "umap", group.by = "orig.ident", label = FALSE, pt.size = 0.1) + 
  theme_classic(base_size = 11) + labs(title = "B. Alignment by Visium Slide") + NoLegend()

p3_umap <- DimPlot(combined_seurat, reduction = "umap", group.by = "Race", pt.size = 0.1, cols = c("Black" = "#E64B35FF", "White" = "#4DBBD5FF")) + 
  theme_classic(base_size = 11) + labs(title = "C. Alignment by Cohort Race")

layout_umap_fixed <- (p1_umap | p2_umap | p3_umap) + 
  plot_annotation(title = "Cohort Global Integration Summary", theme = theme(plot.title = element_text(face="bold", size=14, hjust=0.5)))

ggsave(filename = paste0(out_dir, "P4_Global_UMAP_Panel.pdf"), plot = layout_umap_fixed, width = 16, height = 5)

# Save the physical spatial plots to disk safely using a PDF driver
pdf(paste0(out_dir, "P5_Spatial_Dim_.pdf"), width = 20, height = 30)
SpatialDimPlot(combined_seurat, ncol = 4, pt.size.factor = 1.2)
dev.off()

pdf(paste0(out_dir, "P5_Spatial_feature.pdf"), width = 30, height = 30)
SpatialFeaturePlot(combined_seurat, features = c("BRCA1", "BRCA2"), ncol = 4, pt.size.factor = 1.2)
dev.off()


# ------------------------------------------------------------------------------
# STEP 8: Downstream Global Marker Identification & Race Comparisons
# ------------------------------------------------------------------------------
combined_seurat <- JoinLayers(combined_seurat, assay = "SCT")

niche_markers <- FindAllMarkers(
  object = combined_seurat, 
  assay = "SCT", 
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25,
  verbose = TRUE
)
write.csv(niche_markers, paste0(out_dir, "Cohort_Global_Niche_Markers.csv"), row.names = FALSE)

# Generate demographic tracking identities
combined_seurat$cluster_race <- paste0("Cluster", combined_seurat$seurat_clusters, "_", combined_seurat$Race)
Idents(combined_seurat) <- "cluster_race"

# Compare groups within Niche Cluster 0 using corrected covariates parameter
cluster0_race_diff <- FindMarkers(
  object = combined_seurat,
  ident.1 = "Cluster0_Black",
  ident.2 = "Cluster0_White",
  assay = "SCT",
  latent.vars = "sample_id")

write.csv(cluster0_race_diff, paste0(out_dir, "Black_vs_White_Differential_Expression.csv"))
# ==============================================================================



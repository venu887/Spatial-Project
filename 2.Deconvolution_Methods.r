# ==============================================================================
# PIPELINE: RCTD CELL TYPE DECONVOLUTION ON VISIUM SPATIAL DATA
# ==============================================================================
rm(list=ls())
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(corrplot)
library(GEOquery)
library(spacexr)

out_dir <- "/home/mekalav/Spatial/output/"
dir_path <- "/home/mekalav/Spatial/GSE176078"
ext_path <- file.path(dir_path, "extracted/Wu_etal_2021_BRCA_scRNASeq")

# Set global multi-threading options for compilation/download safety
options(timeout = 600)
Sys.setenv(MAKEFLAGS = "-j2")

# ------------------------------------------------------------------------------
# 1. DOWNLOAD, EXTRACT, AND PREPARE SINGLE-CELL REFERENCE (Wu et al. GSE176078)
# ------------------------------------------------------------------------------
message(">>> Loading Wu et al. Single-Cell Reference Data...")

if (!file.exists(file.path(dir_path, "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"))) {
  getGEOSuppFiles("GSE176078", fetch_files = TRUE, filter_regex = "BRCA_scRNASeq", baseDir = dir_path)
}

if (!dir.exists(ext_path)) {
  tar_file <- file.path(dir_path, "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz")
  untar(tar_file, exdir = file.path(dir_path, "extracted"))
}

counts <- ReadMtx(
  mtx = file.path(ext_path, "count_matrix_sparse.mtx"),
  cells = file.path(ext_path, "count_matrix_barcodes.tsv"),
  features = file.path(ext_path, "count_matrix_genes.tsv"),
  feature.column = 1
)

metadata <- read.csv(file.path(ext_path, "metadata.csv"), row.names = 1)

# Create reference Seurat object
sc_ref <- CreateSeuratObject(counts = counts, meta.data = metadata)

# Extract major cell types
cell_types <- setNames(as.factor(sc_ref$celltype_major), colnames(sc_ref))

# Filter rare cell populations (< 25 cells) to prevent RCTD errors
cell_counts <- table(cell_types)
valid_types <- names(cell_counts[cell_counts >= 25])
keep_indices <- cell_types %in% valid_types

# Apply filtering
filtered_cell_types <- factor(cell_types[keep_indices])
raw_counts <- GetAssayData(sc_ref, assay = "RNA", layer = "counts")[, keep_indices]
nUMI <- colSums(raw_counts)

# Build and save the RCTD Reference Object
message(">>> Building RCTD single-cell reference framework...")
rctd_reference <- Reference(counts = raw_counts, cell_types = filtered_cell_types, nUMI = nUMI)

saveRDS(rctd_reference, file.path(dir_path, "matched_breast_scRNA_reference.rds"))
# Free memory
rm(sc_ref, counts, raw_counts, metadata); gc()

# ------------------------------------------------------------------------------
# 2. LOAD VISIUM SPATIAL DATA AND BUILD SPATIALRNA OBJECT (USING SCT ASSAY)
# -----------------------------------------------------------------------------
combined_seurat <- readRDS(file.path(out_dir, "integrated_data1.rds"))

# Join layers 
combined_seurat <- JoinLayers(combined_seurat, assay = "SCT")
combined_seurat <- JoinLayers(combined_seurat, assay = "Spatial")

# Extract raw counts
visium_counts <- GetAssayData(combined_seurat, assay = "SCT", layer = "counts")

# Extract and reformat tissue coordinates (must have columns 'x' and 'y')
visium_coords <- GetTissueCoordinates(combined_seurat)
spatial_coords <- data.frame(
  x = visium_coords$x,
  y = visium_coords$y,
  row.names = rownames(visium_coords))

spatial_nUMI <- colSums(visium_counts)

# Build SpatialRNA Puck Object
puck <- SpatialRNA(
  coords = spatial_coords,
  counts = visium_counts,
  nUMI   = spatial_nUMI
)

# ------------------------------------------------------------------------------
# 3. INITIALIZE AND EXECUTE RCTD DECONVOLUTION
# ------------------------------------------------------------------------------
my_RCTD <- create.RCTD(
  spatialRNA = puck,
  reference  = rctd_reference,
  max_cores  = 4)

# Run RCTD in full mode for Visium multi-cell spots
my_RCTD <- run.RCTD(my_RCTD, doublet_mode = "full")

# ------------------------------------------------------------------------------
# 4. IMPORT PROPORTIONS BACK INTO SEURAT OBJECT AND SAVE
# ------------------------------------------------------------------------------
# Extract weights and normalize rows to proportions summing to 1
weights_matrix <- as.matrix(my_RCTD@results$weights)
norm_weights <- sweep(weights_matrix, 1, rowSums(weights_matrix), "/")
# Create new Assay containing cell-type proportions (Seurat v5 matrix format)
proportions_assay <- CreateAssayObject(counts = t(norm_weights))
# Inject into Seurat object FIRST, then set as DefaultAssay
combined_seurat[["RCTD"]] <- proportions_assay
DefaultAssay(combined_seurat) <- "RCTD"

# Save updated integrated dataset
saveRDS(combined_seurat, file.path(out_dir, "final_integrated_deconvolved.rds"))
message(">>> Saved successfully to final_integrated_deconvolved.rds")

# ------------------------------------------------------------------------------
# 5. GENERATE SPATIAL FEATURE VISUALIZATIONS
# ------------------------------------------------------------------------------
# Check available deconvolved cell types
available_celltypes <- rownames(combined_seurat[["RCTD"]])
print(available_celltypes)

# Define cell types to map
target_profiles <- c("Myeloid", "T-cells")

pdf(file.path(out_dir, "Fig6_Spatial_Cell_Type_Deconvolution_Maps.pdf"), width = 14, height = 12)
SpatialFeaturePlot(
  object = combined_seurat,
  features = target_profiles,
  images = c("B01", "W01"), # Direct visual comparison between Black and White cohort samples
  ncol = 2,
  pt.size.factor = 1.5,
  alpha = 0.9
) + plot_annotation(title = "Comparative Spatial Cell Type Infiltration In Situ")
dev.off()






# ==============================================================================
# CARD (Spatially-Informed Deconvolution)
# CARD utilizes spatial (x,y) coordinates to account for spatial dependencies between neighboring Visium spots.
# ==============================================================================
library(CARD)
library(Seurat)

# 1. Extract raw spatial counts and tissue coordinates
spatial_counts <- GetAssayData(combined_seurat, assay = "Spatial", layer = "counts")
visium_coords <- GetTissueCoordinates(combined_seurat)
spatial_coords <- data.frame(
  x = visium_coords$x,
  y = visium_coords$y,
  row.names = rownames(visium_coords)
)

# 2. Extract single-cell reference counts and metadata
sc_counts <- GetAssayData(sc_ref_seurat, assay = "RNA", layer = "counts")
sc_meta <- data.frame(
  cellID = colnames(sc_ref_seurat),
  cellType = sc_ref_seurat$cell_type,
  sampleID = "sc_reference",
  row.names = colnames(sc_ref_seurat)
)

# 3. Create CARD Object and run deconvolution
CARD_obj <- createCARDObject(
  sc_count = sc_counts,
  sc_meta = sc_meta,
  spatial_count = spatial_counts,
  spatial_location = spatial_coords,
  ct.varname = "cellType",
  ct.select = NULL,
  sample.varname = "sampleID"
)

CARD_obj <- CARD_deconv(CARD_obj)

# 4. Inject cell type proportions into Seurat
card_weights <- CARD_obj@Proportion_ND
combined_seurat[["CARD"]] <- CreateAssayObject(counts = t(card_weights))






# ==============================================================================
# SPOTlight (NMF + Non-Negative Least Squares)
# SPOTlight identifies marker genes for each cell type in your single-cell reference and uses Non-Negative Matrix Factorization (NMF) to deconvolute spatial spots.
# ==============================================================================
library(SPOTlight)
library(Seurat)

# 1. Find marker genes for single-cell reference cell types
Idents(sc_ref_seurat) <- "cell_type"
sc_markers <- FindAllMarkers(
  sc_ref_seurat, 
  only.pos = TRUE, 
  logfc.threshold = 0.25, 
  min.pct = 0.25
)

# 2. Extract spatial counts matrix
spatial_counts <- GetAssayData(combined_seurat, assay = "Spatial", layer = "counts")

# 3. Run SPOTlight deconvolution
spotlight_res <- spotlight_deconv(
  se_sc = sc_ref_seurat,
  counts_spatial = spatial_counts,
  mks = sc_markers,
  cluster_markers = "cluster",
  cl_type = "cell_type"
)

# 4. Extract proportion matrix and inject into Seurat
spotlight_weights <- spotlight_res$mat
# Set spot barcodes as row names if missing
rownames(spotlight_weights) <- colnames(combined_seurat)

combined_seurat[["SPOTlight"]] <- CreateAssayObject(counts = t(spotlight_weights))









# ==============================================================================
# SpatialDecon (Log-Normal Regression)
# SpatialDecon uses log-normal regression to fit single-cell reference expression profiles to spatial counts, making it resistant to background noise.
# ==============================================================================
library(SpatialDecon)
library(Seurat)

# 1. Build average gene expression matrix per cell type from reference
sc_counts <- GetAssayData(sc_ref_seurat, assay = "RNA", layer = "counts")
cell_types <- as.factor(sc_ref_seurat$cell_type)

# Create mean expression matrix (genes x cell_types)
ref_matrix <- make.neyt.matrix(
  counts = as.matrix(sc_counts),
  celltypes = cell_types
)

# 2. Extract raw spatial count matrix
spatial_counts <- as.matrix(GetAssayData(combined_seurat, assay = "Spatial", layer = "counts"))

# 3. Execute SpatialDecon
spatialdecon_res <- spatialdecon(
  norm = spatial_counts,
  bg = 1, # Background noise cutoff
  X = ref_matrix
)

# 4. Inject cell type fractions into Seurat
spatialdecon_weights <- t(spatialdecon_res$prop_of_all)
combined_seurat[["SpatialDecon"]] <- CreateAssayObject(counts = spatialdecon_weights)




# ==============================================================================
# Seurat Native Anchor Transfer (FindTransferAnchors)
# Seurat's built-in label transfer uses Canonical Correlation Analysis (CCA) anchors to map single-cell identity probabilities onto spatial spots.
# ==============================================================================
library(Seurat)

# 1. Ensure reference is normalized with SCTransform
sc_ref_seurat <- SCTransform(sc_ref_seurat, verbose = FALSE)

# 2. Find anchors between reference and spatial object
anchors <- FindTransferAnchors(
  reference = sc_ref_seurat,
  query = combined_seurat,
  normalization.method = "SCT",
  recompute.residuals = FALSE
)

# 3. Transfer cell type probability scores
predictions <- TransferData(
  anchorset = anchors,
  refdata = sc_ref_seurat$cell_type,
  prediction.assay.name = "SeuratTransfer",
  weight.reduction = combined_seurat[["pca"]],
  dims = 1:30
)

# 4. Add predictions directly to Seurat object
combined_seurat[["SeuratTransfer"]] <- predictions




#@@@@@@@@@@@@@@@@@@@@@@@@@ Comparing Visualizations Across Methods
# Compare T-cell predictions across RCTD, CARD, and SPOTlight
p1 <- SpatialFeaturePlot(combined_seurat, features = "RCTD_T-cells", images = "B01") + ggtitle("RCTD")
p2 <- SpatialFeaturePlot(combined_seurat, features = "CARD_T-cells", images = "B01") + ggtitle("CARD")
p3 <- SpatialFeaturePlot(combined_seurat, features = "SPOTlight_T-cells", images = "B01") + ggtitle("SPOTlight")

p1 | p2 | p3






















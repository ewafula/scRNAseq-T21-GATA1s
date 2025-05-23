# Integrate samples
 
# Eric Wafula
# 2024

# Load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(glmGamPoi))
suppressPackageStartupMessages(library(future))

# Activate parallelization
# plan("multiprocess", workers = 30)
options(future.globals.maxSize= 400000000000)

# Functions
create_seurat_object <- function(sample_data_dir, sample_name) {
  # create seurat object
  sample_data <- Seurat::Read10X(data.dir = sample_data_dir)
  seurat_object <- Seurat::CreateSeuratObject(counts = sample_data, 
                                              project = sample_name, 
                                              min.cells = 3, 
                                              min.features = 200)
  return(seurat_object)
}


# set seed for reproducibility
set.seed(123)

# establish directories
data_dir <- file.path("data")
plots_dir <- file.path("plots", "primary")
results_dir <- file.path("results", "primary")
RDS_dir <- file.path("seurat_objects")

# Create results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Create RDS folder if it doesn't exist
if (!dir.exists(seurat_objects)) {
  dir.create(seurat_objects)
}

# create seurat objects from individual samples
T21.wtGATA1.iPSC <- create_seurat_object(file.path(data_dir, "time_course", "DS145_EM9_D7", 
                                                   "doubletfinder"), "T21wtGATA1iPSC")
T21.GATA1s.iPSC <- create_seurat_object(file.path(data_dir, "time_course", "TMD145_SCM71_D7", 
                                                  "doubletfinder"), "T21GATA1siPSC")
Euploid.wtGATA1.iPSC <- create_seurat_object(file.path(data_dir, "time_course", "TMD145_EM331434_D7", 
                                                       "doubletfinder"), "EuploidwtGATA1iPSC")
Euploid.GATA1s.iPSC <- create_seurat_object(file.path(data_dir, "time_course", "TMD145_EM3319_D7", 
                                                      "doubletfinder"), "EuploidGATA1siPSC")
T21.wtGATA1.FL21 <- create_seurat_object(file.path(data_dir, "primary", "HTS_SO48_01_FL21", 
                                                   "doubletfinder"), "T21wtGATA1FL21")
T21.wtGATA1.FL33 <- create_seurat_object(file.path(data_dir, "primary", "HTS_SO48_02_FL33", 
                                                   "doubletfinder"), "T21wtGATA1FL33")
T21.PennFL32 <- create_seurat_object(file.path(data_dir, "primary", "HTS_SO92_13_FL32_T21_CD34plus", 
                                               "doubletfinder"), "T21PennFL32")
Euploid.wtGATA1.FL38 <- create_seurat_object(file.path(data_dir, "primary", "HTS_SO48_03_FL38", 
                                                       "doubletfinder"), "EuploidwtGATA1FL38")
Euploid.wtGATA1.FL50 <- create_seurat_object(file.path(data_dir, "primary", "HTS_SO48_04_FL50", 
                                                       "doubletfinder"), "EuploidwtGATA1FL50")
Euploid.wtGATA1.FL54 <- create_seurat_object(file.path(data_dir, "primary", "HTS_SO48_05_FL54", 
                                                       "doubletfinder"), "EuploidwtGATA1FL54")
T21.GATA1s.TMD145 <- create_seurat_object(file.path(data_dir, "primary", "HTS_SO48_06_TMD145", 
                                                    "doubletfinder"), "T21GATA1sTMD145")
T21.GATA1s.TMD160 <- create_seurat_object(file.path(data_dir, "primary", "HTS_SO48_07_TMD160", 
                                                    "doubletfinder"), "T21GATA1sTMD160")

# merge seurat sample objects
merged_obj <- merge(T21.wtGATA1.iPSC, y = c(T21.GATA1s.iPSC, Euploid.wtGATA1.iPSC, 
                                            Euploid.GATA1s.iPSC, T21.wtGATA1.FL21,
                                            T21.wtGATA1.FL33, T21.PennFL32, Euploid.wtGATA1.FL38,
                                            Euploid.wtGATA1.FL50, Euploid.wtGATA1.FL54,
                                            T21.GATA1s.TMD145, T21.GATA1s.TMD160), 
                    add.cell.ids = c("T21wtGATA1iPSC", "T21GATA1siPSC", "EuploidwtGATA1iPSC",
                                     "EuploidGATA1siPSC", "T21wtGATA1FL21", "T21wtGATA1FL33",
                                     "T21PennFL32", "EuploidwtGATA1FL38", "EuploidwtGATA1FL50",
                                     "EuploidwtGATA1FL54", "T21GATA1sTMD145", "T21GATA1sTMD160"),
                               project = "GATA1")

# remove objects to conserve memory
rm(T21.wtGATA1.iPSC, T21.GATA1s.iPSC, Euploid.wtGATA1.iPSC, Euploid.GATA1s.iPSC,
   T21.wtGATA1.FL21, T21.wtGATA1.FL33, T21.PennFL32, Euploid.wtGATA1.FL38, 
   Euploid.wtGATA1.FL50, Euploid.wtGATA1.FL54, T21.GATA1s.TMD145, T21.GATA1s.TMD160)

# format metadata
merged_obj$sample <- rownames(merged_obj@meta.data)
merged_obj@meta.data <- merged_obj@meta.data %>% 
  dplyr::mutate(sample = gsub("-1", "", sample)) %>% 
  tidyr::separate(col = sample, into = c("sample", "barcode"), sep = "_") 

# apply sctransform normalization and pca dimension reduction 
merged_obj <- 
  Seurat::SCTransform(merged_obj, variable.features.n = 6000)
merged_obj <- Seurat::RunPCA(merged_obj, npcs = 50)
pdf(file=file.path(plots_dir, "timecourse-d7-and-primary-elbowplot.pdf"))
print(Seurat::ElbowPlot(merged_obj, ndims = 50))
dev.off()

# integrate samples with previously used integration method for SCTransform 
# in Seurat v4 integration
merged_obj <- 
  Seurat::IntegrateLayers(object = merged_obj,
                          method = RPCAIntegration, orig.reduction = "pca",
                          normalization.method = "SCT", verbose = FALSE,
                          new.reduction = "integrated.rpca")
integrated_obj <- 
  Seurat::FindNeighbors(merged_obj, dims = 1:20, reduction = "integrated.rpca")
integrated_obj <-  Seurat::FindClusters(integrated_obj, resolution = 0.4)
integrated_obj <- Seurat::RunUMAP(integrated_obj, dims = 1:20, 
                                  reduction = "integrated.rpca")

# remove object to conserve memory
rm(merged_obj)

# Save seurat object to RDS file
file_name <- file.path(RDS_dir, "timecourse-d7-and-primary.RDS")
saveRDS(integrated_obj, file = file_name)

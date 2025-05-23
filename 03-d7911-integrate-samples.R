# Integrate samples
 
# Eric Wafula
# 10/20/2024

# Load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(glmGamPoi))
suppressPackageStartupMessages(library(DoubletFinder))
suppressPackageStartupMessages(library(DropletUtils))
suppressPackageStartupMessages(library(future))

# Activate parallelization
# plan("multicore", workers = 10)
options(future.globals.maxSize = 200000*1024^2, future.seed = NULL)


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

# establish base, data, results, and plots directories
root_dir <- here::here()
data_dir <- file.path(root_dir, "data", "timecourse")
plots_dir <- file.path(root_dir, "plots", "timecourse")
RDS_dir <- file.path(root_dir, "seurat_objects")
if (!dir.exists(RDS_dir)) {
  dir.create(RDS_dir)
}

# assign variables 
# percent_mito = 10
ncps = 20
nfeatures = 6000

# create seurat objects from individual samples
T21.wtGATA1.D7 <- create_seurat_object(file.path(data_dir, "DS145_EM9_D7", 
                                          "doubletfinder"), "T21wtGATA1D7")
T21.GATA1s.D7 <- create_seurat_object(file.path(data_dir, "TMD145_SCM71_D7", 
                                              "doubletfinder"), "T21GATA1sD7")
Euploid.wtGATA1.D7 <- create_seurat_object(file.path(data_dir, "TMD145_EM331434_D7", 
                                              "doubletfinder"), "EuploidwtGATA1D7")
Euploid.GATA1s.D7 <- create_seurat_object(file.path(data_dir, "TMD145_EM3319_D7", 
                                             "doubletfinder"), "EuploidGATA1sD7")
T21.wtGATA1.D9 <- create_seurat_object(file.path(data_dir, "DS145_EM9_D9", 
                                                 "doubletfinder"), "T21wtGATA1D9")
T21.GATA1s.D9 <- create_seurat_object(file.path(data_dir, "TMD145_SCM71_D9", 
                                                "doubletfinder"), "T21GATA1sD9")
Euploid.wtGATA1.D9 <- create_seurat_object(file.path(data_dir, "TMD145_EM331434_D9", 
                                                     "doubletfinder"), "EuploidwtGATA1D9")
Euploid.GATA1s.D9 <- create_seurat_object(file.path(data_dir, "TMD145_EM3319_D9", 
                                                    "doubletfinder"), "EuploidGATA1sD9")
T21.wtGATA1.D11 <- create_seurat_object(file.path(data_dir, "DS145_EM9_D11", 
                                                 "doubletfinder"), "T21wtGATA1D11")
T21.GATA1s.D11 <- create_seurat_object(file.path(data_dir, "TMD145_SCM71_D11", 
                                                "doubletfinder"), "T21GATA1sD11")
Euploid.wtGATA1.D11 <- create_seurat_object(file.path(data_dir, "TMD145_EM331434_D11", 
                                                     "doubletfinder"), "EuploidwtGATA1D11")
Euploid.GATA1s.D11 <- create_seurat_object(file.path(data_dir, "TMD145_EM3319_D11", 
                                                    "doubletfinder"), "EuploidGATA1sD11")

# merge seurat sample objects
merged_seurat_objects <- merge(T21.wtGATA1.D7, 
                               y = c(T21.GATA1s.D7, Euploid.wtGATA1.D7, Euploid.GATA1s.D7, 
                                     T21.wtGATA1.D9, T21.GATA1s.D9, Euploid.wtGATA1.D9, 
                                     Euploid.GATA1s.D9, T21.wtGATA1.D11, T21.GATA1s.D11, 
                                     Euploid.wtGATA1.D11, Euploid.GATA1s.D11), 
                               add.cell.ids = c("T21wtGATA1D7", "T21GATA1sD7", 
                                                "EuploidwtGATA1D7", "EuploidGATA1sD7",
                                                "T21wtGATA1D9", "T21GATA1sD9",
                                                "EuploidwtGATA1D9", "EuploidGATA1sD9",
                                                "T21wtGATA1D11", "T21GATA1sD11",
                                                "EuploidwtGATA1D11", "EuploidGATA1sD11"),
                               project = "GATA1")

# remove objects to conserve memory
rm(T21.wtGATA1.D7, T21.GATA1s.D7, Euploid.wtGATA1.D7, Euploid.GATA1s.D7,
   T21.wtGATA1.D9, T21.GATA1s.D9, Euploid.wtGATA1.D9, Euploid.GATA1s.D9,
   T21.wtGATA1.D11, T21.GATA1s.D11, Euploid.wtGATA1.D11, Euploid.GATA1s.D11)

# format metadata
merged_seurat_objects$Sample <- rownames(merged_seurat_objects@meta.data)
merged_seurat_objects@meta.data <- merged_seurat_objects@meta.data %>% 
  tidyr::separate(col = Sample, into = c("Sample", "Barcode"))

# filter mitochondria
# not needed for doubletfinder or soupx filtered datasets
# samples filtered already filtered to remove mitochondria
merged_seurat_objects[["percent.mt"]] <- 
  Seurat::PercentageFeatureSet(merged_seurat_objects, pattern = "^MT[-\\.]")
# merged_seurat_objects <- subset(merged_seurat_objects, percent.mt < percent_mito)

# apply sctransform normalization
seurat_objects.list <- Seurat::SplitObject(merged_seurat_objects, 
                                           split.by = "Sample")

# remove object to conserve memory
rm(merged_seurat_objects)

for (i in 1:length(seurat_objects.list)) { 
  seurat_objects.list[[i]] <- 
    Seurat::SCTransform(seurat_objects.list[[i]], method = "glmGamPoi", 
                        variable.features.n = nrow(seurat_objects.list[[i]]), 
                        ncells = ncol(seurat_objects.list[[i]]))
  seurat_objects.list[[i]] <- Seurat::RunPCA(seurat_objects.list[[i]], 
                                             npcs = ncps,)
  seurat_objects.list[[i]] <- Seurat::FindNeighbors(seurat_objects.list[[i]],
                                                    dims = 1:ncps)
  seurat_objects.list[[i]] <- Seurat::FindClusters(seurat_objects.list[[i]], 
                                                   resolution = 0.4)
  seurat_objects.list[[i]] <- Seurat::RunUMAP(seurat_objects.list[[i]],
                                              dims = 1:ncps)
}

# reference samples to void the following error when running FindIntegrationAnchors()
# Error in .M2C(newTMat(i = c(ij1[, 1], ij2[, 1]), j = c(ij1[, 2], ij2[,  :
# unable to aggregate TsparseMatrix with 'i' and 'j' slots of length exceeding 2^31-1
# Calls: <Anonymous> ... FindIntegrationMatrix -> - -> - -> .Arith.Csparse -> .M2C
# Execution halted
reference_dataset <- which(names(seurat_objects.list) %in% 
                             c("EuploidwtGATA1D7", "EuploidwtGATA1D9", "EuploidwtGATA1D11"))

# integrate samples using seurat rpca
features <- Seurat::SelectIntegrationFeatures(object.list = seurat_objects.list, 
                                              nfeatures = nfeatures)
seurat_objects.list <- Seurat::PrepSCTIntegration(object.list = seurat_objects.list, 
                                                  anchor.features = features)
sample.anchors <- Seurat::FindIntegrationAnchors(object.list = seurat_objects.list, 
                                                 normalization.method = "SCT",
                                                 anchor.features = features, 
                                                 dims = 1:ncps, reduction = "rpca", 
                                                 k.anchor = 20, reference = reference_dataset)

# remove object to conserve memory
rm(seurat_objects.list)

# all_genes <- lapply(seurat_objects.list, row.names) %>% Reduce(intersect, .)
samples.integrated.sct <- Seurat::IntegrateData(anchorset = sample.anchors, 
                                                normalization.method = "SCT", 
                                                #features.to.integrate = all_genes) # memory failure
                                                features.to.integrate = features)

# remove object to conserve memory
rm(sample.anchors)

# run a single analysis on integrated samples
samples.integrated.sct <- Seurat::RunPCA(samples.integrated.sct, ncps = 50,
                                         verbose = FALSE)
pdf(file=file.path(plots_dir, "timecourse-d7911-elbowplot.pdf"))
print(Seurat::ElbowPlot(samples.integrated.sct, ndims = 50, reduction = "pca"))
dev.off()
samples.integrated.sct <- Seurat::RunUMAP(samples.integrated.sct, reduction = "pca",
                                          dims = 1:ncps)
samples.integrated.sct <- Seurat::FindNeighbors(samples.integrated.sct, 
                                                reduction = "pca",
                                                dims = 1:ncps)
samples.integrated.sct <- Seurat::FindClusters(samples.integrated.sct, 
                                               resolution = 
                                                 c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
saveRDS(samples.integrated.sct, file = file.path(RDS_dir, "timecourse-d7911.RDS"))

# remove object to conserve memory
rm(samples.integrated.sct)

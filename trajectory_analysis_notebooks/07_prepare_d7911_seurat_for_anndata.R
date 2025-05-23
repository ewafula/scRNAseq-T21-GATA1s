# Prepare data from D7/D9/D11 timecourse Seurat object for converting
# to Python anndata

# Eric Wafula
# 2024

# Load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(Matrix))


# set seed for reproducibility
set.seed(123)

# establish directories
root_dir <- here::here()
data_dir <- file.path(root_dir, "..", "data", "timecourse")
RDS_dir <- file.path(root_dir, "seurat_objects")
results_dir <- file.path(root_dir, "results", "timecourse")

# load integrated and  annotated Seurat RDS Object
integrated <- readRDS(file.path(RDS_dir, "timecourse-d7911-renamed.RDS"))

# save metadata table:
integrated$barcode <- colnames(integrated)
integrated_metadata <- integrated@meta.data %>% 
  dplyr::select(-Sample, -Barcode, -integrated_snn_res.0.1, 
                -integrated_snn_res.0.2, -integrated_snn_res.0.3, 
                -integrated_snn_res.0.5, -integrated_snn_res.0.6, 
                -integrated_snn_res.0.7, -nCount_SCT, -nFeature_SCT,
                -SCT_snn_res.0.4) %>% 
  # tidyr::separate(col = barcode, into = c("sample", "barcode"), sep = "_") %>%
  # dplyr::select(-sample) %>% 
  dplyr::mutate(condition = case_when(grepl("T21", orig.ident) ~ "T21", 
                                      grepl("Euploid", orig.ident) ~ "Euploid"), 
                genotype = case_when(grepl("wtGATA1", orig.ident) ~ "wtGATA1", 
                                     grepl("GATA1s", orig.ident) ~ "GATA1s"), 
                time = case_when(grepl("D7", orig.ident) ~ "D7", 
                                 grepl("D9", orig.ident) ~ "D9", 
                                 grepl("D11", orig.ident) ~ "D11")) %>% 
  readr::write_csv(file.path(results_dir,
                             "timecourse-d7911-renamed-metadata.csv"))

# write integrated UMAP embeddings ti file
integrated@reductions$umap@cell.embeddings %>% as.data.frame() %>% 
  rownames_to_column(var = "barcode") %>% 
  readr::write_csv(file.path(results_dir,
                             "timecourse-d7911-integrated-umap-embeddings.csv"))

# join sample count Layers
integrated <- SeuratObject::JoinLayers(integrated, assay = "RNA")

# write expression counts matrix
integrated_count_matrix <- Seurat::GetAssayData(integrated, assay="RNA", layer = "counts")
Matrix::writeMM(d7911_count_matrix, 
                file=file.path(results_dir,
                               "timecourse-d7911-renamed-counts.mtx"))

# write dimensionality reduction pca matrix
integrated_pca_matrix <- integrated@reductions$pca@cell.embeddings %>%
  as.data.frame() %>% 
  readr::write_csv(file.path(results_dir,
                             "timecourse-d7911-renamed-pca.csv"))

# write gene names
integrated_count_matrix %>% rownames() %>% 
  readr::write_lines(file.path(results_dir,
                               "timecourse-d7911-renamed-genes.txt"))


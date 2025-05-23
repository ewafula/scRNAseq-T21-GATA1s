# Remove doublets and multiplets cell with additional filtering of mitochondria
# and features doublets are identified in the 10X data and the doublet barcodes 
# are removed from the soupX data. 

# Eric Wafula
# 10/20/2024


# Load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(DoubletFinder))
suppressPackageStartupMessages(library(DropletUtils))

# functions
create_seurat_object <- function(sample_data_dir, sample_name) {
  # create seurat object
  sample_data <- Seurat::Read10X(data.dir = sample_data_dir)
  seurat_object <- Seurat::CreateSeuratObject(counts = sample_data, 
                                              project = sample_name, 
                                              min.cells = 3, 
                                              min.features = 200)
  return(seurat_object)
}


filter_seurat_object <- function(seurat_object, components, percent_mito) {
  # create seurat object, filter mitochondria, normalize and cluster
  seurat_object[["percent.mt"]] <- 
    Seurat::PercentageFeatureSet(seurat_object, pattern = "^MT[-\\.]")
  seurat_object <- subset(seurat_object, percent.mt < percent_mito)
  seurat_object <- Seurat::NormalizeData(seurat_object)
  seurat_object <- Seurat::FindVariableFeatures(seurat_object, 
                                                selection.method = "vst",
                                                nfeatures = 2000)
  seurat_object <- Seurat::ScaleData(seurat_object)
  seurat_object <- Seurat::RunPCA(seurat_object)
  seurat_object <- Seurat::FindNeighbors(seurat_object, dims = 1:components)
  seurat_object <- Seurat::FindClusters(seurat_object)
  seurat_object <- Seurat::RunUMAP(seurat_object, dims = 1:components)
  return(seurat_object)
}


get_pK <- function(seurat_object, sample_name, results_dir) {
  # calculate pk
  sweep.res.list <- DoubletFinder::paramSweep_v3(seurat_object, 
                                                 PCs = 1:components, 
                                                 sct = FALSE)
  sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- DoubletFinder::find.pK(sweep.stats)
  # write pk values to file
  file_name <- here::here(results_dir, paste0(sample_name,'_pk_values.txt'))
  write.table(bcmvn, file_name, append = FALSE, sep = '\t', dec = '.', 
              quote=FALSE, row.names = FALSE, col.names = TRUE)
  # sweep plot
  pdf(file=here::here(plots_dir, paste0(sample_name, '_pk_sweep_plot.pdf')))
  pK=as.numeric(as.character(bcmvn$pK))
  BCmetric=bcmvn$BCmetric
  pK_choose = pK[which(BCmetric %in% max(BCmetric))]
  par(mar=c(5,4,4,1)+.1,cex.main=1.2,font.main=2)
  plot(x = pK, y = BCmetric, pch = 16, type="b",
       col = "blue",lty=1, cex = 1.2)
  abline(v=pK_choose,lwd=2,col='red',lty=2)
  title("The BCmvn distributions")
  text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 2,col = "red")
  dev.off()
  # get the pK values 
  pK <- readr::read_tsv(file_name) %>% 
    dplyr::filter(BCmetric == max(BCmetric)) %>% 
    dplyr::pull(pK)
  return(pK)
}

run_doubletfinder <- function(seurat_object, doublet_rate, pN, pK) {
  pK <- as.numeric(pK)
  annotations <- seurat_object@meta.data$seurat_clusters
  homotypic.prop <- DoubletFinder::modelHomotypic(annotations)
  nExp_poi <- round(doublet_rate*nrow(seurat_object@meta.data))  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  seurat_object <- 
    DoubletFinder::doubletFinder_v3(seurat_object, PCs = 1:components, pN = pN, 
                                    pK = pK, nExp = nExp_poi.adj, 
                                    reuse.pANN = FALSE, sct = FALSE)
  # umap plot
  pdf(file=here::here(plots_dir, paste0(sample_name, '_umap_doublets.pdf')))
  metadata <- names(seurat_object@meta.data)
  metadata <- metadata[grep("^DF\\.classifications_", metadata)]
  print(Seurat::DimPlot(seurat_object, reduction = "umap", group.by = metadata))
  dev.off()
  return(seurat_object)
}

remove_doublets <- function(seurat_object, doubletfinder_dir, percent_mito, 
                            sample_name, soupx_dir) {
  # get doublets
  metadata <- names(seurat_object@meta.data)
  metadata <- metadata[grep("^DF\\.classifications_", metadata)]
  doublets <- seurat_object[, seurat_object@meta.data[, metadata] == "Doublet"]
  doublet_ids <- unique(colnames(doublets))
  # Create seurat object from soupX data
  seurat_object_soupx <- create_seurat_object(soupx_dir, sample_name)
  # Filter mitochondria
  seurat_object_soupx[["percent.mt"]] <- 
    Seurat::PercentageFeatureSet(seurat_object_soupx, pattern = "^MT[-\\.]")
  seurat_object_soupx <- subset(seurat_object_soupx, percent.mt < percent_mito)
  # Remove doublets from soupX Seurat object
  seurat_object_soupx <- 
    seurat_object_soupx[,!colnames(seurat_object_soupx) %in% doublet_ids]
  # write doublet filtered results
  DropletUtils::write10xCounts(x = seurat_object_soupx@assays$RNA@counts, 
                               path = doubletfinder_dir, overwrite = TRUE, version='3')
}


# set up optparse options
option_list <- list(
  make_option(opt_str = "--sample_dir", type = "character", default = NULL,
              help = "Directory with 10x cellranger filtered or SoupX samples data",
              metavar = "character"),
  make_option(opt_str = "--sample_name", type = "character", default = NULL,
              help = "Sample name - sub-directory with sample cellranger results",
              metavar = "character"),
  make_option(opt_str = "--percent_mito", type = "numeric", default = NULL,
              help = "Maximum percent ratio of mitochondria cell content to allow",
              metavar = "character"),
  make_option(opt_str = "--components", type = "numeric", default = NULL,
              help = "Number of components for PCA",
              metavar = "character"),
  make_option(opt_str = "--doublet_rate", type = "numeric", default = NULL,
              help = "Doublet rate from recovered cells report in the assay kit manual",
              metavar = "character"),
  make_option(opt_str = "--soupx", action = "store_true", default = TRUE,
              help = "Remove doublets from soupX filtered data",
              metavar = "character")
)

# parse parameter options
opt <- parse_args(OptionParser(option_list = option_list))
sample_dir <- opt$sample_dir
sample_name <- opt$sample_name
percent_mito <- opt$percent_mito
components <- opt$components
doublet_rate <- opt$doublet_rate
soupx <- opt$soupx

# assign variables
pN <- 0.25 # default recommended for DoubletFinder

# set seed for reproducibility
set.seed(123)

# establish base, data, results, and plots directories
root_dir <- here::here()
data_dir <- file.path(root_dir, "data")
plots_dir <- file.path(root_dir, "plots", "doubletfinder")
results_dir <- file.path(root_dir, "results", "doubletfinder")
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# set input and output sub-directories
cellranger_dir <- file.path(sample_dir, "outs", "filtered_feature_bc_matrix")
soupx_dir <- file.path(sample_dir, "soupx")
doubletfinder_dir <- file.path(sample_dir, "doubletfinder")
if (!dir.exists(doubletfinder_dir)) {
  dir.create(doubletfinder_dir)
}

# run doublet removal analysis 
seurat_object <- create_seurat_object(cellranger_dir, sample_name)
seurat_object_filtered <- filter_seurat_object(seurat_object, components, 
                                               percent_mito)
pK <- get_pK(seurat_object_filtered, sample_name, results_dir)
seurat_object_doublets <- run_doubletfinder(seurat_object_filtered, 
                                            doublet_rate, pN, pK)
if (soupx) {
  remove_doublets(seurat_object_doublets, doubletfinder_dir, percent_mito,
                  sample_name, soupx_dir)
} else {
  # Remove doublets from filtered Seurat object
  metadata <- names(seurat_object_doublets@meta.data)
  metadata <- metadata[grep("^DF\\.classifications_", metadata)]
  doublets <- 
    seurat_object_doublets[, seurat_object_doubletst@meta.data[, metadata] == "Doublet"]
  doublet_ids <- unique(colnames(doublets))
  seurat_object_doublets <- 
    seurat_object_doublets[,!colnames(seurat_object_doublets) %in% doublet_ids]
  # write doublet filtered results
  DropletUtils::write10xCounts(x = seurat_object_doublets@assays$RNA@counts, 
                               path = doubletfinder_dir, overwrite = TRUE, version='3')
}
# Removed pK sweep plot created by doubletfinder.A better pK sweep plot plot 
# already created using the pK distribution values output file.
unlink("Rplots.pdf")

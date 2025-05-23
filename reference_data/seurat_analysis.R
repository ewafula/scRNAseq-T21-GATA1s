library(Seurat)
library(dplyr)
library(Matrix)
require("ggrepel")

# Load the data
load("Rdata/seurat_object.Rdata")

# update load old seurat object to accommodate new features
sample <- Seurat::UpdateSeuratObject(sample)

# Count Mito genes and make plots
mito.features = grep(pattern="^MT-", x=rownames(x=sample), value=TRUE)
percent.mito = Matrix::colSums(x=GetAssayData(object=sample, slot="counts")[mito.features,]) /
  Matrix::colSums(x=GetAssayData(object=sample, slot="counts"))
sample[["percent.mito"]] = percent.mito
VlnPlot(object=sample, features=c("nFeature_RNA","nCount_RNA","percent.mito"),
        group.by="orig.ident", ncol=3, pt.size=0.1)
ggsave("plots/quality_by_origin.pdf", width=12, height=6.8)
FeatureScatter(object=sample, feature1="nCount_RNA", feature2="percent.mito")
FeatureScatter(object=sample, feature1="nCount_RNA", feature2="nFeature_RNA")

# Filter cells and Normalize data
sample = subset(x=sample, subset=nFeature_RNA>100 & percent.mito<0.1)
sample = NormalizeData(object = sample, normalization.method = "LogNormalize", 
                                   scale.factor = 10000)

# Find variable Genes and scale data by number of UMIs and Mito gene percentage
sample = FindVariableFeatures(object=sample, selection.method="vst")
cc_genes = c(cc.genes$s.genes, cc.genes$g2m.genes)
sample@assays$RNA@var.features = setdiff(sample@assays$RNA@var.features, cc_genes)
length(sample@assays$RNA@var.features)
sample = ScaleData(object=sample)

# Run PCA
sample = RunPCA(object=sample, verbose=FALSE)
DimPlot(sample, reduction="pca", group.by="ident")
ElbowPlot(object=sample, ndims=30)

# Run TSNE
pca_dims = 1:20
sample = RunTSNE(object=sample, dims=pca_dims, verbose=FALSE)
sample = RunUMAP(object=sample, dims=pca_dims, return.model = TRUE, verbose=FALSE)
sample = FindNeighbors(object=sample, dims=pca_dims, verbose=FALSE)
sample = FindClusters(object=sample, resolution=0.5, verbose=FALSE)

# Make DM plots
set.seed(1)
random_order = sample(colnames(sample))
DimPlot(sample, label=TRUE, reduction="umap", group.by="ident", pt.size=0.1)
ggsave("plots/umap_clusterRX105.pdf", width=7, height=6)
DimPlot(sample, label=FALSE, reduction="umap", group.by="orig.ident", pt.size=0.1, cells=random_order)
ggsave("plots/umap_originRX105.pdf", width=9.5, height=6)
DimPlot(sample, label=TRUE, reduction="tsne", group.by="ident", pt.size=0.1)
ggsave("plots/tsne_clusterRX105.pdf", width=7, height=6)
DimPlot(sample, label=FALSE, reduction="tsne", group.by="orig.ident", pt.size=0.1, cells=random_order)
ggsave("plots/tsne_originRX105.pdf", width=9.5, height=6)

# Save the Seurat object
save(file="Rdata/seurat_object_rerun.Rdata", sample)

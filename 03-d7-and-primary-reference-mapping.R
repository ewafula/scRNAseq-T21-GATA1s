# Reference based mapping 

# Eric Wafula
# 09/09/2024

# Load libraries
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
require("ggrepel")

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

set.seed(123)

# Load the data
load("reference_data/Rdata/seurat_object.Rdata")

# update load old seurat object to accommodate new features
reference <- Seurat::UpdateSeuratObject(sample)

# Rename identity classes
reference <-
  Seurat::RenameIdents(reference, "0" = "Monocyte/Macrophage", "1" = "Lympho-Myeloid",
                       "2" = "Erythroid", "3" = "Megakaryocytic", "4" = "Granulocyte",
                       "5" = "Erythroid", "6" = "Monocyte/Macrophage", "7" = "Lymphoid",
                       "8" = "Erythroid", "9" = "HSCs", "10" = "Monocyte/Macrophage",
                       "11" = "Monocyte/Macrophage", "12" = "Megakaryocytic",
                       "13" = "Granulocyte", "14" = "HSCs", "15" = "Monocyte/Macrophage",
                       "16" = "B Lymphocyte", "17" = "Lymphoid", "18" = "Monocyte/Macrophage",
                       "19" = "HSCs", "20" = "Monocyte/Macrophage", "21" = "T Lymphocyte",
                       "22" = "B Lymphocyte")


Idents(reference) <-  
  factor(Idents(reference),
         levels = c("Monocyte/Macrophage", "Lympho-Myeloid", "Erythroid", "Megakaryocytic",
                    "Granulocyte", "Lymphoid", "HSCs", "B Lymphocyte", "T Lymphocyte"))

reference[["cell_type"]] <- Seurat::Idents(reference)

# plot UMAP grouped by renamed cell types (idents)
Seurat::DimPlot(reference, reduction = "umap", label=FALSE, pt.size=0.1)
ggsave("plots/reference/reference_umap_cell_types.pdf", width=9.5, height=6)

# plot UMAP grouped by reference sample
Seurat::DimPlot(reference, reduction = "umap", label=FALSE, pt.size=0.1, 
                group.by = "orig.ident")
ggsave("plots/reference/reference_umap_samples.pdf", width=9.5, height=6)

### transfer cell type annotation to rerun Seurat v5 object metadata 

# get annotated cell types
cell_types <- reference@meta.data %>% dplyr::select(cell_type) %>% 
  tibble::rownames_to_column(var = "barcode")

# Load the rerun data
rm(sample)
load("reference_data/Rdata/seurat_object_rerun.Rdata")

# add annotated cell types
barcodes = rownames(sample@meta.data)
cell_types = cell_types[match(barcodes, cell_types$barcode), ] #match the order of seurat barcode
sample@meta.data$cell_type = cell_types$cell_type #put cell types into the metadata
reference <- sample 

# save reference seurat object
saveRDS(reference, file = file.path("seurat_objects/reference.RDS"))

####  prepare and load query samples

# query sample dataset directory
data_dir <- "data"

# load query samples 
T21.wtGATA1.D7 <- 
  create_seurat_object(file.path(data_dir, "timecourse", "DS145_EM9_D7", 
                                 "doubletfinder"), "T21wtGATA1D7")
T21.GATA1s.D7 <- 
  create_seurat_object(file.path(data_dir, "timecourse", "TMD145_SCM71_D7", 
                                 "doubletfinder"), "T21GATA1sD7")
Euploid.wtGATA1.D7 <- 
  create_seurat_object(file.path(data_dir, "timecourse", "TMD145_EM331434_D7", 
                                 "doubletfinder"), "EuploidwtGATA1D7")
Euploid.GATA1s.D7 <- 
  create_seurat_object(file.path(data_dir, "timecourse", "TMD145_EM3319_D7", 
                                 "doubletfinder"), "EuploidGATA1sD7")
T21.wtGATA1.FL21 <- 
  create_seurat_object(file.path(data_dir, "primary", "HTS_SO48_01_FL21", 
                                 "doubletfinder"), "T21wtGATA1FL21")
T21.wtGATA1.FL33 <- 
  create_seurat_object(file.path(data_dir, "primary", "HTS_SO48_02_FL33", 
                                 "doubletfinder"), "T21wtGATA1FL33")
T21.PennFL32 <- 
  create_seurat_object(file.path(data_dir, "primary", "HTS_SO92_13_FL32_T21_CD34plus", 
                                 "doubletfinder"), "T21PennFL32")
Euploid.wtGATA1.FL38 <- 
  create_seurat_object(file.path(data_dir, "primary", "HTS_SO48_03_FL38", 
                                 "doubletfinder"), "EuploidwtGATA1FL38")
Euploid.wtGATA1.FL50 <- 
  create_seurat_object(file.path(data_dir, "primary", "HTS_SO48_04_FL50", 
                                 "doubletfinder"), "EuploidwtGATA1FL50")
Euploid.wtGATA1.FL54 <- 
  create_seurat_object(file.path(data_dir, "primary", "HTS_SO48_05_FL54",
                                 "doubletfinder"), "EuploidwtGATA1FL54")
T21.GATA1s.TMD145 <-
  create_seurat_object(file.path(data_dir, "primary", "HTS_SO48_06_TMD145",
                                 "doubletfinder"), "T21GATA1sTMD145")
T21.GATA1s.TMD160 <- 
  create_seurat_object(file.path(data_dir, "primary", "HTS_SO48_07_TMD160", 
                                 "doubletfinder"), "T21GATA1sTMD160")

# merge seurat sample objects
query <- merge(T21.wtGATA1.D7, y = c(T21.GATA1s.D7, Euploid.wtGATA1.D7, Euploid.GATA1s.D7,
                                     T21.wtGATA1.FL21, T21.wtGATA1.FL33, T21.PennFL32,
                                     Euploid.wtGATA1.FL38, Euploid.wtGATA1.FL50, Euploid.wtGATA1.FL54, 
                                     T21.GATA1s.TMD145, T21.GATA1s.TMD160), 
               add.cell.ids = c("T21wtGATA1D7", "T21GATA1sD7", "EuploidwtGATA1D7", "EuploidGATA1sD7", 
                                "T21wtGATA1FL21", "T21wtGATA1FL33", "T21PennFL32",
                                "EuploidwtGATA1FL38", "EuploidwtGATA1FL50", "EuploidwtGATA1FL54", 
                                "T21GATA1sTMD145", "T21GATA1sTMD160"), project = "GATA1")

# join layers
query <- SeuratObject::JoinLayers(query)

# remove objects to conserve memory
rm(T21.wtGATA1.D7, T21.GATA1s.D7, Euploid.wtGATA1.D7, Euploid.GATA1s.D7,
   T21.wtGATA1.FL21, T21.wtGATA1.FL33, T21.PennFL32,
   Euploid.wtGATA1.FL38, Euploid.wtGATA1.FL50, Euploid.wtGATA1.FL54, 
   T21.GATA1s.TMD145, T21.GATA1s.TMD160)

# format metadata
query$sample <- rownames(query@meta.data)
query@meta.data <- query@meta.data %>% 
  tidyr::separate(col = sample, into = c("sample", "barcode"), sep = "_")

# Filter cells and Normalize data
query = NormalizeData(object = query, normalization.method = "LogNormalize", 
                       scale.factor = 10000)

# get reference anchors
reference_anchors <- Seurat::FindTransferAnchors(reference = reference, 
                                                 query = query, dims = 1:30,
                                                 reference.reduction = "pca")

# query to reference mapping
query <- Seurat::MapQuery(anchorset = reference_anchors, reference = reference,
                           query = query, refdata = list(cell_type = "cell_type"), 
                           reference.reduction = "pca",reduction.model = "umap")

# reference and query sampple UMAPs
cols = c('Monocyte/Macrophage'='#F68282', 'Lympho-Myeloid'='#31C53F', 'HSCs'='#8494FF',
         'Erythroid'='#B95FBB', 'Megakaryocytic'='#ff9a36','Granulocyte'='#28CECA',
         'Lymphoid'='#AC8F14', 'B Lymphocyte'='#25aff5', 'T Lymphocyte'='#E6C122')

p1 <- DimPlot(reference, reduction = "umap", group.by = "cell_type", label = TRUE, 
              label.size = 3, repel = TRUE, cols = cols) + NoLegend() + 
  ggtitle("Reference annotations")
p2 <- DimPlot(query, reduction = "ref.umap", group.by = "predicted.cell_type", 
              label = TRUE, label.size = 3, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("Query transferred labels")
p1 + p2
ggsave("plots/reference/reference_query_umap_samples_primary_plus_d7.pdf", width=12, height=6)

# add tissue type metadata to query
query@meta.data <- query@meta.data %>% 
  dplyr::mutate(sample_type = 
                  case_when(orig.ident %in% 
                              c("T21wtGATA1D7", "T21GATA1sD7", "EuploidwtGATA1D7",
                                "EuploidGATA1sD7") ~ "iPSC",
                            orig.ident %in% 
                              c("T21wtGATA1FL21", "T21wtGATA1FL33",  "T21PennFL32", 
                                "EuploidwtGATA1FL38", "EuploidwtGATA1FL50", 
                                "EuploidwtGATA1FL54") ~ "Fetal Liver",
                            orig.ident %in% 
                              c("T21GATA1sTMD145", "T21GATA1sTMD160") ~ "PBMC" )) %>% 
  dplyr::mutate(sample = 
                  case_when(sample == "T21wtGATA1D7" ~ "T21/wtGATA1/iPSC/D7", 
                            sample == "T21GATA1sD7" ~ "T21/GATA1s/iPSC/D7", 
                            sample == "EuploidwtGATA1D7" ~ "Euploid/wtGATA1/iPSC/D7", 
                            sample == "EuploidGATA1sD7" ~ "Euploid/GATA1s/iPSC/D7", 
                            sample == "T21wtGATA1FL21" ~ "T21/wtGATA1/FL/W22",
                            sample == "T21wtGATA1FL33" ~ "T21/wtGATA1/FL/W15",
                            sample == "T21PennFL32" ~ "T21/FL/W16",
                            sample == "EuploidwtGATA1FL38" ~ "Euploid/wtGATA1/FL/W15",
                            sample == "EuploidwtGATA1FL50" ~ "Euploid/wtGATA1/FL/W15/D3",
                            sample == "EuploidwtGATA1FL54" ~ "Euploid/wtGATA1/FL/W22",
                            sample == "T21GATA1sTMD145" ~ "T21/GATA1s/TMD/68%",
                            sample == "T21GATA1sTMD160" ~ "T21/GATA1s/TMD/22%"))

# save query seurat object
saveRDS(query, file = file.path("seurat_objects/query_primary_plus_d7.RDS"))

# query sample types UMAPs
p1 <- DimPlot(query, reduction = "ref.umap", group.by = "predicted.cell_type", 
              label = TRUE, label.size = 3, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("Query transferred labels")
p2 <- DimPlot(query, reduction = "ref.umap", group.by = "sample_type", 
               label = FALSE, label.size = 3) + 
  ggtitle("Sample types")
p1 + p2
ggsave("plots/reference/reference_query_umap_sample_types_primary_plus_d7.pdf", 
       width=15, height=6)
p2
ggsave("plots/reference/query_umap_sample_types_primary_plus_d7.pdf",
       width=9.5, height=6)

# reference split UMAPs by samples 
AGM_4wk_658 <- subset(reference, subset = orig.ident == "AGM_4wk_658")
p1 <- DimPlot(AGM_4wk_658, reduction = "umap", group.by = "cell_type", 
              label = FALSE, label.size = 5, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("AGM/CS14/4wk") +
  theme(plot.title = element_text(size = 40, face = "bold"))
AGM_5wk_555 <- subset(reference, subset = orig.ident == "AGM_5wk_555")
p2 <- DimPlot(AGM_5wk_555, reduction = "umap", group.by = "cell_type", 
              label = FALSE, label.size = 5, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("AGM/CS15a/5wk") +
  theme(plot.title = element_text(size = 40, face = "bold"))
AGM_5wk_575 <- subset(reference, subset = orig.ident == "AGM_5wk_575")
p3 <- DimPlot(AGM_5wk_575, reduction = "umap", group.by = "cell_type", 
              label = FALSE, label.size = 5, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("AGM/CS15b/5wk") +
  theme(plot.title = element_text(size = 40, face = "bold"))
AGM_6wk_563 <- subset(reference, subset = orig.ident == "AGM_6wk_563")
p4 <- DimPlot(AGM_6wk_563, reduction = "umap", group.by = "cell_type", 
              label = FALSE, label.size = 5, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("AGM/CS17/6wk") +
  theme(plot.title = element_text(size = 40, face = "bold"))
Liver_4wk_658 <- subset(reference, subset = orig.ident == "Liver_4wk_658")
p5 <- DimPlot(Liver_4wk_658, reduction = "umap", group.by = "cell_type", 
              label = FALSE, label.size = 5, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("Liver/CS14/4wk") +
  theme(plot.title = element_text(size = 40, face = "bold"))
Liver_5wk_575 <- subset(reference, subset = orig.ident == "Liver_5wk_575")
p6 <- DimPlot(Liver_5wk_575, reduction = "umap", group.by = "cell_type", 
              label = FALSE, label.size = 5, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("Liver/CS15b/5wk") +
  theme(plot.title = element_text(size = 40, face = "bold"))
Liver_6wk_563 <- subset(reference, subset = orig.ident == "Liver_6wk_563")
p7 <- DimPlot(Liver_6wk_563, reduction = "umap", group.by = "cell_type", 
              label = FALSE, label.size = 5, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("Liver/CS17/6wk") +
  theme(plot.title = element_text(size = 40, face = "bold"))
Liver_8wk_553 <- subset(reference, subset = orig.ident == "Liver_8wk_553")
p8 <- DimPlot(Liver_8wk_553, reduction = "umap", group.by = "cell_type", 
              label = FALSE, label.size = 5, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("Liver/8wk") +
  theme(plot.title = element_text(size = 40, face = "bold"))
Liver_11wk_569 <- subset(reference, subset = orig.ident == "Liver_11wk_569")
p9 <- DimPlot(Liver_11wk_569, reduction = "umap", group.by = "cell_type", 
              label = FALSE, label.size = 5, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("Liver/11wk") +
  theme(plot.title = element_text(size = 40, face = "bold"))
Liver_15wk_101 <- subset(reference, subset = orig.ident == "Liver_15wk_101")
p10 <- DimPlot(Liver_15wk_101, reduction = "umap", group.by = "cell_type", 
               label = FALSE, label.size = 5, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("Liver/15wk") +
  theme(plot.title = element_text(size = 40, face = "bold"))
design <- "ABC#
           DEF#
           GHIJ"
p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 +
  patchwork::plot_layout(design = design)
ggsave("plots/reference/reference_umap_samples_split.pdf", width=46, height=36)

# replot UMAP grouped by reference sample
DimPlot(reference, reduction = "umap", group.by = "cell_type", 
        label = FALSE, label.size = 3, cols = cols) + ggtitle("All Samples")
ggsave("plots/reference/reference_umap_cell_types_v2.pdf", width=9.5, height=9.5)

# query split UMPAPs by samples 
T21wtGATA1D7 <- subset(query, subset = sample == "T21/wtGATA1/iPSC/D7")
p1 <- DimPlot(T21wtGATA1D7, reduction = "ref.umap", group.by = "predicted.cell_type", 
              raster = FALSE, label = FALSE, label.size = 3, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("T21/wtGATA1/iPSC/D7") + 
  theme(plot.title = element_text(size = 25, face = "bold"))
T21GATA1sD7 <- subset(query, subset = sample == "T21/GATA1s/iPSC/D7")
p2 <- DimPlot(T21GATA1sD7, reduction = "ref.umap", group.by = "predicted.cell_type", 
              raster = FALSE, label = FALSE, label.size = 3, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("T21/GATA1s/iPSC/D7") +
  theme(plot.title = element_text(size = 25, face = "bold"))
EuploidwtGATA1D7 <- subset(query, subset = sample == "Euploid/wtGATA1/iPSC/D7")
p3 <- DimPlot(EuploidwtGATA1D7, reduction = "ref.umap", group.by = "predicted.cell_type", 
              raster = FALSE, label = FALSE, label.size = 3, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("Euploid/wtGATA1/iPSC/D7") +
  theme(plot.title = element_text(size = 25, face = "bold"))
EuploidGATA1sD7 <- subset(query, subset = sample == "Euploid/GATA1s/iPSC/D7")
p4 <- DimPlot(EuploidGATA1sD7, reduction = "ref.umap", group.by = "predicted.cell_type", 
              raster = FALSE, label = FALSE, label.size = 3, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("Euploid/GATA1s/iPSC/D7") +
  theme(plot.title = element_text(size = 25, face = "bold"))
T21wtGATA1FL15 <- subset(query, subset = sample == "T21/wtGATA1/FL/W15")
p5 <- DimPlot(T21wtGATA1FL15, reduction = "ref.umap", group.by = "predicted.cell_type", 
              raster = FALSE, label = FALSE, label.size = 3, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("T21/wtGATA1/FL/W15") +
  theme(plot.title = element_text(size = 25, face = "bold"))
T21PennFL32 <- subset(query, subset = sample == "T21/FL/W16")
p6 <- DimPlot(T21PennFL32, reduction = "ref.umap", group.by = "predicted.cell_type", 
              raster = FALSE, label = FALSE, label.size = 3, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("T21/FL/W16") +
  theme(plot.title = element_text(size = 25, face = "bold"))
T21wtGATA1FL22 <- subset(query, subset = sample == "T21/wtGATA1/FL/W22")
p7 <- DimPlot(T21wtGATA1FL22, reduction = "ref.umap", group.by = "predicted.cell_type", 
              raster = FALSE, label = FALSE, label.size = 3, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("T21/wtGATA1/FL/W22") +
  theme(plot.title = element_text(size = 25, face = "bold"))
EuploidwtGATA1FL15 <- subset(query, subset = sample == "Euploid/wtGATA1/FL/W15")
p8 <- DimPlot(EuploidwtGATA1FL15, reduction = "ref.umap", group.by = "predicted.cell_type", 
              raster = FALSE, label = FALSE, label.size = 3, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("Euploid/wtGATA1/FL/W15") +
  theme(plot.title = element_text(size = 25, face = "bold"))
EuploidwtGATA1FL153 <- subset(query, subset = sample == "Euploid/wtGATA1/FL/W15/D3")
p9 <- DimPlot(EuploidwtGATA1FL153, reduction = "ref.umap", group.by = "predicted.cell_type", 
              raster = FALSE, label = FALSE, label.size = 3, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("Euploid/wtGATA1/FL/W15/D3") +
  theme(plot.title = element_text(size = 25, face = "bold"))
EuploidwtGATA1FL22 <- subset(query, subset = sample == "Euploid/wtGATA1/FL/W22")
p10 <- DimPlot(EuploidwtGATA1FL22, reduction = "ref.umap", group.by = "predicted.cell_type", 
               raster = FALSE, label = FALSE, label.size = 3, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("Euploid/wtGATA1/FL/W22") +
  theme(plot.title = element_text(size = 25, face = "bold"))
T21GATA1sTMD22 <- subset(query, subset = sample == "T21/GATA1s/TMD/22%")
p11 <- DimPlot(T21GATA1sTMD22, reduction = "ref.umap", group.by = "predicted.cell_type", 
               raster = FALSE, label = FALSE, label.size = 3, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("T21/GATA1s/TMD/22%") +
  theme(plot.title = element_text(size = 25, face = "bold"))
T21GATA1sTMD68 <- subset(query, subset = sample == "T21/GATA1s/TMD/68%")
p12 <- DimPlot(T21GATA1sTMD68, reduction = "ref.umap", group.by = "predicted.cell_type", 
               raster = FALSE, label = FALSE, label.size = 3, repel = TRUE, cols = cols) + 
  NoLegend() +  ggtitle("T21/GATA1s/TMD/68%") +
  theme(plot.title = element_text(size = 25, face = "bold"))
(p1 | p2 | p5 | p6 | p7 | p11) / (p3 | p4 |  p8 | p9 | p10 | p12)
ggsave("plots/reference/query_umap_samples_split_primary_plus_d7.pdf", 
       width=72, height=26, limitsize = FALSE)

# write sample cell counts per cell type to file
query@meta.data %>%  dplyr::select(sample, predicted.cell_type) %>% 
  dplyr::count(sample, predicted.cell_type) %>%  
  tidyr::pivot_wider(names_from = predicted.cell_type, values_from = n) %>% 
  replace(is.na(.), 0) %>% 
  readr::write_tsv("results/reference/primary_plus_and_d7_query_sample_cell_type_counts.tsv")


# Sample cluster proportion 
tt<-prop.table(table(query$predicted.cell_type, query$sample), margin = 2)
df<- reshape::melt(tt)
colnames(df)<-c("CellType", "Sample", "Proportion")
df$CellType<-factor(df$CellType, levels = rev(c("HSCs", "Erythroid", "Megakaryocytic", "Lympho-Myeloid",
                                                "Granulocyte", "Monocyte/Macrophage", "Lymphoid", "B Lymphocyte", 
                                                "T Lymphocyte")))
df$Sample<-factor(df$Sample, levels = rev(c("Euploid/wtGATA1/iPSC/D7", "Euploid/GATA1s/iPSC/D7",
                                            "T21/wtGATA1/iPSC/D7" ,"T21/GATA1s/iPSC/D7",
                                            "Euploid/wtGATA1/FL/W15", "Euploid/wtGATA1/FL/W15/D3",
                                            "Euploid/wtGATA1/FL/W22", "T21/wtGATA1/FL/W15", "T21/FL/W16", 
                                            "T21/wtGATA1/FL/W22", "T21/GATA1s/TMD/22%", "T21/GATA1s/TMD/68%")))
pdf(file = "plots/reference/query_umap_predicted_cell_types_primary_plus_d7-proportions.pdf")
p<-ggplot(df, aes(Proportion, Sample, fill = CellType)) + geom_bar(stat = "identity") +
  scale_fill_manual(values = cols, name = "Cell Type", guide = guide_legend(reverse = TRUE)) + 
  theme_classic() +
  theme(axis.text = element_text(size = 8), axis.line.x = element_blank(),
        axis.line.y = element_blank(), axis.title.y = element_blank(), 
        aspect.ratio = .8)
print(p)
dev.off()

# plot UMAP grouped by predicted.cell_type
DimPlot(query, reduction = "ref.umap", group.by = "predicted.cell_type",
              label = FALSE, label.size = 3, cols = cols) + ggtitle("All Samples")
ggsave("plots/reference/query_umap_predicted_cell_types_primary_plus_d7.pdf", width=9.5, height=6)


# gene sets scorcards 
scorecard1 <- c("MEIS2", "GATA1", "GATA2", "RUNX1", "GYPA", "GPI", "LMO2", 
                "SPI1", "SPN", "MYC", "MYB")
scorecard1 <- rev(factor(scorecard1))

scorecard2 <- c("MEIS2", "GATA1", "GATA2", "RUNX1", "ETV6", "ITGA2B", "GYPA", 
                "HBG1", "KLF1", "LMO2", "SPI1", "MYC", "MYB")
scorecard2 <- rev(factor(scorecard2))

#### sample ordering scorecard1
# Erythroid expression
erys <- subset(query, subset = predicted.cell_type == "Erythroid")
erys@meta.data$sample <- factor(erys@meta.data$sample,
                                levels = c("Euploid/wtGATA1/iPSC/D7", "Euploid/GATA1s/iPSC/D7",
                                           "T21/wtGATA1/iPSC/D7" ,"T21/GATA1s/iPSC/D7",
                                           "Euploid/wtGATA1/FL/W15", "Euploid/wtGATA1/FL/W15/D3",
                                           "Euploid/wtGATA1/FL/W22", "T21/wtGATA1/FL/W15", "T21/FL/W16", 
                                           "T21/wtGATA1/FL/W22", "T21/GATA1s/TMD/22%", "T21/GATA1s/TMD/68%"))
DotPlot(erys, features=scorecard1, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/reference/query_umap_samples_primary_plus_d7_scorecard1_ery.pdf", height = 6)

# Megakaryocyte expression
mks <- subset(query, subset = predicted.cell_type == "Megakaryocytic")
# sample ordering scorecard1
mks@meta.data$sample <- factor(mks@meta.data$sample,
                               levels = c("Euploid/wtGATA1/iPSC/D7", "Euploid/GATA1s/iPSC/D7",
                                          "T21/wtGATA1/iPSC/D7" ,"T21/GATA1s/iPSC/D7",
                                          "Euploid/wtGATA1/FL/W15", "Euploid/wtGATA1/FL/W15/D3",
                                          "Euploid/wtGATA1/FL/W22", "T21/wtGATA1/FL/W15", "T21/FL/W16", 
                                          "T21/wtGATA1/FL/W22", "T21/GATA1s/TMD/22%", "T21/GATA1s/TMD/68%"))
DotPlot(mks , features=scorecard1, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/reference/query_umap_samples_primary_plus_d7_scorecard1_mk.pdf", height = 6)


#### sample ordering scorecard2
# Erythroid expression
erys@meta.data$sample <- factor(erys@meta.data$sample,
                                levels = c("Euploid/wtGATA1/iPSC/D7","Euploid/wtGATA1/FL/W15",
                                           "Euploid/wtGATA1/FL/W15/D3", "Euploid/wtGATA1/FL/W22", 
                                           "Euploid/GATA1s/iPSC/D7", "T21/wtGATA1/iPSC/D7", 
                                           "T21/wtGATA1/FL/W15", "T21/FL/W16", "T21/wtGATA1/FL/W22",
                                           "T21/GATA1s/iPSC/D7", "T21/GATA1s/TMD/22%", "T21/GATA1s/TMD/68%"))
DotPlot(erys, features=scorecard2, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/reference/query_umap_samples_primary_plus_d7_scorecard2_ery.pdf", height = 7)

# Megakaryocyte expression
mks@meta.data$sample <- factor(mks@meta.data$sample,
                               levels = c("Euploid/wtGATA1/iPSC/D7","Euploid/wtGATA1/FL/W15",
                                          "Euploid/wtGATA1/FL/W15/D3", "Euploid/wtGATA1/FL/W22", 
                                          "Euploid/GATA1s/iPSC/D7", "T21/wtGATA1/iPSC/D7", 
                                          "T21/wtGATA1/FL/W15", "T21/FL/W16", "T21/wtGATA1/FL/W22",
                                          "T21/GATA1s/iPSC/D7", "T21/GATA1s/TMD/22%", "T21/GATA1s/TMD/68%"))
DotPlot(mks, features=scorecard2, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/reference/query_umap_samples_primary_plus_d7_scorecard2_mk.pdf", height = 7)


####### combine scorecards 
scorecard <- unique(c(scorecard1, scorecard2))

#### sample ordering1
# Erythroid expression
erys <- subset(query, subset = predicted.cell_type == "Erythroid")
erys@meta.data$sample <- factor(erys@meta.data$sample,
                                levels = c("Euploid/wtGATA1/iPSC/D7", "Euploid/GATA1s/iPSC/D7",
                                           "T21/wtGATA1/iPSC/D7" ,"T21/GATA1s/iPSC/D7",
                                           "Euploid/wtGATA1/FL/W15", "Euploid/wtGATA1/FL/W15/D3",
                                           "Euploid/wtGATA1/FL/W22", "T21/wtGATA1/FL/W15", "T21/FL/W16", 
                                           "T21/wtGATA1/FL/W22", "T21/GATA1s/TMD/22%", "T21/GATA1s/TMD/68%"))
DotPlot(erys, features=scorecard, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/reference/query_umap_samples_primary_plus_d7-ordering1-ery.pdf", height = 7)

# Megakaryocyte expression
mks <- subset(query, subset = predicted.cell_type == "Megakaryocytic")

mks@meta.data$sample <- factor(mks@meta.data$sample,
                               levels = c("Euploid/wtGATA1/iPSC/D7", "Euploid/GATA1s/iPSC/D7",
                                          "T21/wtGATA1/iPSC/D7" ,"T21/GATA1s/iPSC/D7",
                                          "Euploid/wtGATA1/FL/W15", "Euploid/wtGATA1/FL/W15/D3",
                                          "Euploid/wtGATA1/FL/W22", "T21/wtGATA1/FL/W15", "T21/FL/W16", 
                                          "T21/wtGATA1/FL/W22", "T21/GATA1s/TMD/22%", "T21/GATA1s/TMD/68%"))
DotPlot(mks , features=scorecard, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/reference/query_umap_samples_primary_plus_d7-ordering1-mk.pdf", height = 7)


#### sample ordering2
# Erythroid expression
erys@meta.data$sample <- factor(erys@meta.data$sample,
                                levels = c("Euploid/wtGATA1/iPSC/D7","Euploid/wtGATA1/FL/W15",
                                           "Euploid/wtGATA1/FL/W15/D3", "Euploid/wtGATA1/FL/W22", 
                                           "Euploid/GATA1s/iPSC/D7", "T21/wtGATA1/iPSC/D7", 
                                           "T21/wtGATA1/FL/W15", "T21/FL/W16", "T21/wtGATA1/FL/W22",
                                           "T21/GATA1s/iPSC/D7", "T21/GATA1s/TMD/22%", "T21/GATA1s/TMD/68%"))
DotPlot(erys, features=scorecard, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/reference/query_umap_samples_primary_plus_d7-ordering2-ery.pdf", height = 7)

# Megakaryocyte expression
mks@meta.data$sample <- factor(mks@meta.data$sample,
                               levels = c("Euploid/wtGATA1/iPSC/D7","Euploid/wtGATA1/FL/W15",
                                          "Euploid/wtGATA1/FL/W15/D3", "Euploid/wtGATA1/FL/W22", 
                                          "Euploid/GATA1s/iPSC/D7", "T21/wtGATA1/iPSC/D7", 
                                          "T21/wtGATA1/FL/W15", "T21/FL/W16", "T21/wtGATA1/FL/W22",
                                          "T21/GATA1s/iPSC/D7", "T21/GATA1s/TMD/22%", "T21/GATA1s/TMD/68%"))
DotPlot(mks, features=scorecard, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/reference/query_umap_samples_primary_plus_d7-ordering2-mk.pdf", height = 7)

library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(data.table)

# load Seurat RDS Object
res <- readRDS("seurat_objects/timecourse-d7-and-primary.RDS")

# plot clusteringf UMAP
pdf(file = "plots/pub/timecourse-d7-and-primary-umap.pdf")
Seurat::DimPlot(res, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
dev.off()

# plot clusteringf UMAP
pdf(file = "plots/pub/timecourse-d7-and-primary-umap_bigfont.pdf")
plot <- Seurat::DimPlot(res, reduction = "umap") + NoLegend()
LabelClusters(plot, id = "ident", position = "nearest", size = 6, fontface = "bold")
dev.off()

# define marker list for cluster annotation
markers <- list()

# eythroid marker genes
# https://maayanlab.cloud/Harmonizome/gene_set/erythrocyte/TISSUES+Curated+Tissue+Protein+Expression+Evidence+Scores
ery_genes <- c("HBA", "HBB", "HBZ", "KLF1", "HBE1", "GYPA", "GYPB","BBS1", "DMTN", "GSR", 
              "LANCL1", "MPP1", "ATP2B4", "ATP2B1", "CYB5A","HDHD1","PLSCR1",
              "ALDH1A1","CAT","CD47", "BLVRB", "MT1E", "CD44","TSTA3", "ACHE",
              "EPB41","EPB42", "GLO1","BMI1", "GATA1", "ADD2", "PCMT1","CYB5R3",
              "ADD1", "PRDX6", "PRDX2", "PGK1", "STOM", "CD58","CTSE","SKP1",
              "HBG2", "HBA1", "HBA2", "HBG1", "SLC25A37", "ANK1", "TFRC", "CD36",
              "SPTA1", "BPGM", "ALAS2")
markers$Ery <- ery_genes

# megakaryocytes marker genes
# https://maayanlab.cloud/Harmonizome/gene_set/megakaryocyte/TISSUES+Experimental+Tissue+Protein+Expression+Evidence+Scores
mk_genes <- c("PF4", "GP9", "SELP", "DIAPH1", "UGGT2", "PSPC1", "VPS13A", "TMC1",  
              "NBEAL2", "SPTB", "SYTL4", "SPATA17", "KIF16B","MMRN1", "ALB", "GPD2", 
              "UBE2O","ZNF559", "CFAP36", "NCF1", "IRF7","WDFY3", "CASS4", "UTRN", 
              "GPI", "SPTBN4", "CLTC", "WNT8B", "NEXN", "UBR4", "HK1", "HSP90B1", 
              "STIP1", "SPTA1", "ZSCAN12", "TAF1B", "DDI1", "FLNA", "CCDC181", "FGA",
              "SPTBN2","PDIA4", "FGB", "CCT3", "YARS", "VWF", "HMHA1", "ANK1", 
              "GTF2H4", "F5", "EHD1", "NBEAL1", "ACLY", "CREB1", "DAAM2", "P4HB",
              "TPP2", "DAAM1", "LTF", "CNGB1", "DYNC1H1", "VCL", "LRRC57","WDR1", 
              "VCP", "COPA", "PYGB", "PYGL", "CREM", "DCTN1", "WDR44", "SMC3", "UNC13D",
              "PDZD3", "ENPP1", "PIEZO2", "MYO18A", "UGGT1", "SLC12A1", "HSP90AA1",
              "LTBP1", "ITGB3", "ITGA2B","THBS1", "PLEK", "MED12L", "PRKAR2B","ARHGAP6",
              "ESAM", "SERPINB1")
markers$MK <- mk_genes 

# myeloid marker genes #######
# from Chris except the EMP3 and VASP from MsigDB
# https://www.gsea-msigdb.org/gsea/msigdb/cards/MA_MYELOID_DIFFERENTIATION_UP
myeloid_genes <- c("ACSL1", "AZU1", "CEBPA", "CSTA", "CXCL8", "CYBB", "CCL3",
                   "CCL4", "CCL4l2", "CLC", "MPO", "MRC1", "IRAK3", "NRP1", "IKZF2",
                   "PDE4D", "MS4A3", "RBM47", "PRTN3", "PRG2", "PRG3", "SAMHD1",
                   "S100A8", "S100A9", "VIM", "PU1", "SP1", "EMP3", "VASP", "CD74", 
                   "SRGN", "CST3", "AOAH", "SLC9A9", "CTSB", "LGMN", "FTL")
markers$Myeloid <- myeloid_genes

# add scores
res_scores <- UCell::AddModuleScore_UCell(res, features = markers)
res_scores@meta.data <- res_scores@meta.data %>% as.data.frame %>%
  dplyr::rename(Ery = Ery_UCell, MK = MK_UCell, Myeloid = Myeloid_UCell)

# plot marker expression umap and dotplots
p1 <-  Seurat::FeaturePlot(res_scores, reduction = "umap", features = "Ery", label = TRUE) + 
  scale_colour_gradientn(colours = brewer.pal(n = 11, name = "YlOrRd"))
p2 <-  Seurat::FeaturePlot(res_scores, reduction = "umap", features = "MK", label = TRUE) + 
  scale_colour_gradientn(colours = brewer.pal(n = 11, name = "YlOrRd"))
p3 <-  Seurat::FeaturePlot(res_scores, reduction = "umap", features = "Myeloid", label = TRUE) + 
  scale_colour_gradientn(colours = brewer.pal(n = 11, name = "YlOrRd"))
p4 <- Seurat::DotPlot(res_scores, features = c("Ery","MK","Myeloid"), dot.scale = 5) & 
  scale_colour_gradientn(colours = brewer.pal(n = 11, name = "YlOrRd"))
(p1 | p2) / (p3 | p4)
ggsave("plots/pub/timecourse-d7-and-primary-cluster-annotation.pdf", 
       width=12, height=10, limitsize = FALSE)

# Rename identity classes
res_renamed <- Seurat::RenameIdents(res_scores, "0" = "HPC - unbiased", "1" = "HPC - Ery 1",
                                    "2" = "HPC - unbiased", "3" = "HPC - Ery 2", "4" = "HPC - Mk 2", 
                                    "5" = "HPC - Myeloid 1", "6" = "HPC - Myeloid 3",
                                    "7" = "HPC - unbiased", "8" = "HPC - Ery 3", "9" = "HPC - Mk 1", 
                                    "10" = "HPC - Ery 1", "11" = "HPC - unbiased", 
                                    "12" = "HPC - Myeloid 1", "13" = "HPC - unbiased", 
                                    "14" = "HPC - unbiased", "15" = "HPC - unbiased", "16" = "HPC - unbiased",
                                    "17" = "HPC - Mk 3", "18" = "HPC - Myeloid 2", "19" = "HPC - unbiased")

Idents(res_renamed) <-  factor(Idents(res_renamed), 
                               levels = c("HPC - unbiased", 
                                          "HPC - Mk 1", "HPC - Mk 2", "HPC - Mk 3", 
                                          "HPC - Ery 1", "HPC - Ery 2", "HPC - Ery 3", 
                                          "HPC - Myeloid 1", "HPC - Myeloid 2", "HPC - Myeloid 3"))

cols = c("HPC - unbiased" = "snow3", 
         "HPC - Mk 1" = "plum3", "HPC - Mk 2" = "mediumpurple1", "HPC - Mk 3" = "purple4",
         "HPC - Ery 1" = "lightpink3", "HPC - Ery 2" = "violetred3", "HPC - Ery 3" = "red4", 
         "HPC - Myeloid 1" = "turquoise", "HPC - Myeloid 2" = "springgreen3", "HPC - Myeloid 3" = "forestgreen")

pdf(file = "plots/pub/timecourse-d7-and-primary-umap-renamed.pdf", width = 10)
Seurat::DimPlot(res_renamed, reduction = "umap", cols = cols)
dev.off()
pdf(file = "plots/pub/timecourse-d7-and-primary-umap-renamed-v2.pdf")
Seurat::DimPlot(res_renamed, reduction = "umap", cols = cols) + NoLegend()
dev.off()

# add cell type annotations as column in meta data
res_renamed[["cell_type"]] <- Seurat::Idents(res_renamed)

# add tissue type metadata and save object to RDS file
res_renamed[["sample_type"]] <- res_renamed$sample
res_renamed@meta.data <- res_renamed@meta.data %>% 
  dplyr::mutate(sample_type = 
                  case_when(sample_type %in% c("T21wtGATA1iPSC", "T21GATA1siPSC", 
                                               "EuploidwtGATA1iPSC", "EuploidGATA1siPSC") ~ "iPSC",
                            sample_type %in% c("T21wtGATA1FL21", "T21wtGATA1FL33", "T21PennFL32", 
                                               "EuploidwtGATA1FL38", "EuploidwtGATA1FL50", 
                                               "EuploidwtGATA1FL54") ~ "Fetal Liver",
                            sample_type %in% c("T21GATA1sTMD145", "T21GATA1sTMD160") ~ "PBMC")) %>% 
  dplyr::mutate(sample = 
                  case_when(sample == "T21wtGATA1iPSC" ~ "T21/wtGATA1/iPSC/D7", 
                            sample == "T21GATA1siPSC" ~ "T21/GATA1s/iPSC/D7", 
                            sample == "EuploidwtGATA1iPSC" ~ "Euploid/wtGATA1/iPSC/D7", 
                            sample == "EuploidGATA1siPSC" ~ "Euploid/GATA1s/iPSC/D7", 
                            sample == "T21wtGATA1FL21" ~ "T21/wtGATA1/FL/W22",
                            sample == "T21wtGATA1FL33" ~ "T21/wtGATA1/FL/W15",
                            sample == "T21PennFL32" ~ "T21/FL/W16",
                            sample == "EuploidwtGATA1FL38" ~ "Euploid/wtGATA1/FL/W15",
                            sample == "EuploidwtGATA1FL50" ~ "Euploid/wtGATA1/FL/W15/D3",
                            sample == "EuploidwtGATA1FL54" ~ "Euploid/wtGATA1/FL/W22",
                            sample == "T21GATA1sTMD145" ~ "T21/GATA1s/TMD/68%",
                            sample == "T21GATA1sTMD160" ~ "T21/GATA1s/TMD/22%"))
  

pdf(file = "plots/pub/timecourse-d7-and-primary-umap-sample_type.pdf")
Seurat::DimPlot(res_renamed, group.by = "sample_type")
dev.off()

# write cell counts to file
table(res_renamed@meta.data$sample, res_renamed@meta.data$cell_type) %>% 
  as.data.frame.matrix() %>% 
  tibble::rownames_to_column(var = "Sample") %>% 
  readr::write_tsv("results/pub/timecourse-d7-and-primary-umap-renamed-cell-counts.tsv.gz")


# add cell type annotations as column in meta data
res_renamed[["cell_type"]] <- Seurat::Idents(res_renamed)

# query split UMAPs by samples 
T21wtGATA1iPSC <- subset(res_renamed, subset = sample == "T21/wtGATA1/iPSC/D7")
p1 <- DimPlot(T21wtGATA1iPSC, reduction = "umap", group.by = "cell_type", 
              label = FALSE, label.size = 3, repel = TRUE, cols = cols, raster = FALSE) + 
  NoLegend() +  ggtitle("T21/wtGATA1/iPSC/D7") + 
  theme(plot.title = element_text(size = 25, face = "bold"))
T21GATA1siPSC <- subset(res_renamed, subset = sample == "T21/GATA1s/iPSC/D7")
p2 <- DimPlot(T21GATA1siPSC, reduction = "umap", group.by = "cell_type", 
              label = FALSE, label.size = 3, repel = TRUE, cols = cols, raster = FALSE) + 
  NoLegend() +  ggtitle("T21/GATA1s/iPSC/D7") +
  theme(plot.title = element_text(size = 25, face = "bold"))
EuploidwtGATA1iPSC <- subset(res_renamed, subset = sample == "Euploid/wtGATA1/iPSC/D7")
p3 <- DimPlot(EuploidwtGATA1iPSC, reduction = "umap", group.by = "cell_type", 
              label = FALSE, label.size = 3, repel = TRUE, cols = cols, raster = FALSE) + 
  NoLegend() +  ggtitle("Euploid/wtGATA1/iPSC/D7") +
  theme(plot.title = element_text(size = 25, face = "bold"))
EuploidGATA1siPSC <- subset(res_renamed, subset = sample == "Euploid/GATA1s/iPSC/D7")
p4 <- DimPlot(EuploidGATA1siPSC, reduction = "umap", group.by = "cell_type", 
              label = FALSE, label.size = 3, repel = TRUE, cols = cols, raster = FALSE) + 
  NoLegend() +  ggtitle("Euploid/GATA1s/iPSC/D7") +
  theme(plot.title = element_text(size = 25, face = "bold"))
T21wtGATA1FL15 <- subset(res_renamed, subset = sample == "T21/wtGATA1/FL/W15")
p5 <- DimPlot(T21wtGATA1FL15, reduction = "umap", group.by = "cell_type", 
              label = FALSE, label.size = 3, repel = TRUE, cols = cols, raster = FALSE) + 
  NoLegend() +  ggtitle("T21/wtGATA1/FL/W15") +
  theme(plot.title = element_text(size = 25, face = "bold"))
T21PennFL32 <- subset(res_renamed, subset = sample == "T21/FL/W16")
p6 <- DimPlot(T21PennFL32, reduction = "umap", group.by = "cell_type", 
               label = FALSE, label.size = 3, repel = TRUE, cols = cols, raster = FALSE) + 
  NoLegend() +  ggtitle("T21/FL/W16") +
  theme(plot.title = element_text(size = 25, face = "bold"))
T21wtGATA1FL22 <- subset(res_renamed, subset = sample == "T21/wtGATA1/FL/W22")
p7 <- DimPlot(T21wtGATA1FL22, reduction = "umap", group.by = "cell_type", 
              label = FALSE, label.size = 3, repel = TRUE, cols = cols, raster = FALSE) + 
  NoLegend() +  ggtitle("T21/wtGATA1/FL/W22") +
  theme(plot.title = element_text(size = 25, face = "bold"))
EuploidwtGATA1FL15 <- subset(res_renamed, subset = sample == "Euploid/wtGATA1/FL/W15")
p8 <- DimPlot(EuploidwtGATA1FL15, reduction = "umap", group.by = "cell_type", 
              label = FALSE, label.size = 3, repel = TRUE, cols = cols, raster = FALSE) + 
  NoLegend() +  ggtitle("Euploid/wtGATA1/FL/W15") +
  theme(plot.title = element_text(size = 25, face = "bold"))
EuploidwtGATA1FL153 <- subset(res_renamed, subset = sample == "Euploid/wtGATA1/FL/W15/D3")
p9 <- DimPlot(EuploidwtGATA1FL153, reduction = "umap", group.by = "cell_type", 
              label = FALSE, label.size = 3, repel = TRUE, cols = cols, raster = FALSE) + 
  NoLegend() +  ggtitle("Euploid/wtGATA1/FL/W15/D3") +
  theme(plot.title = element_text(size = 25, face = "bold"))
EuploidwtGATA1FL22 <- subset(res_renamed, subset = sample == "Euploid/wtGATA1/FL/W22")
p10 <- DimPlot(EuploidwtGATA1FL22, reduction = "umap", group.by = "cell_type", 
              label = FALSE, label.size = 3, repel = TRUE, cols = cols, raster = FALSE) + 
  NoLegend() +  ggtitle("Euploid/wtGATA1/FL/W22") +
  theme(plot.title = element_text(size = 25, face = "bold"))
T21GATA1sTMD22 <- subset(res_renamed, subset = sample == "T21/GATA1s/TMD/22%")
p11 <- DimPlot(T21GATA1sTMD22, reduction = "umap", group.by = "cell_type", 
               label = FALSE, label.size = 3, repel = TRUE, cols = cols, raster = FALSE) + 
  NoLegend() +  ggtitle("T21/GATA1s/TMD/22%") +
  theme(plot.title = element_text(size = 25, face = "bold"))
T21GATA1sTMD68 <- subset(res_renamed, subset = sample == "T21/GATA1s/TMD/68%")
p12 <- DimPlot(T21GATA1sTMD68, reduction = "umap", group.by = "cell_type", 
               label = FALSE, label.size = 3, repel = TRUE, cols = cols, raster = FALSE) + 
  NoLegend() +  ggtitle("T21/GATA1s/TMD/68%") +
  theme(plot.title = element_text(size = 25, face = "bold"))
(p1 | p2 | p5 | p6 | p7 | p11) / (p3 | p4 |  p8 | p9 | p10 | p12)
ggsave("plots/pub/timecourse-d7-and-primary-umap-renamed-split.pdf", 
       width=72, height=26, limitsize = FALSE)

# Save annotated object to RDS file
file_name <- file.path("seurat_objects/timecourse-d7-and-primary-renamed.RDS")
saveRDS(res_renamed, file = file_name)

# Sample cluster proportion 
tt<-prop.table(table(res_renamed$cell_type, res_renamed$sample), margin = 2)
df<-melt(tt)
colnames(df)<-c("CellType", "Sample", "Proportion")
df$CellType<-factor(df$CellType, levels = rev(c("HPC - unbiased", 
                                                "HPC - Mk 1", "HPC - Mk 2", "HPC - Mk 3", 
                                                "HPC - Ery 1", "HPC - Ery 2", "HPC - Ery 3", 
                                                "HPC - Myeloid 1", "HPC - Myeloid 2", "HPC - Myeloid 3")))
df$Sample<-factor(df$Sample, levels = rev(c("Euploid/wtGATA1/iPSC/D7", "Euploid/GATA1s/iPSC/D7",
                                            "T21/wtGATA1/iPSC/D7" ,"T21/GATA1s/iPSC/D7",
                                            "Euploid/wtGATA1/FL/W15", "Euploid/wtGATA1/FL/W15/D3",
                                            "Euploid/wtGATA1/FL/W22", "T21/wtGATA1/FL/W15", "T21/FL/W16", 
                                            "T21/wtGATA1/FL/W22", "T21/GATA1s/TMD/22%", "T21/GATA1s/TMD/68%")))
pdf(file = "plots/pub/timecourse-d7-and-primary-umap-renamed-proportions.pdf")
p<-ggplot(df, aes(Proportion, Sample, fill = CellType)) + geom_bar(stat = "identity") +
  scale_fill_manual(values = cols, name = "Cell Type", guide = guide_legend(reverse = TRUE)) + 
  theme_classic() +
  theme(axis.text = element_text(size = 8), axis.line.x = element_blank(),
        axis.line.y = element_blank(), axis.title.y = element_blank(), 
        aspect.ratio = .8)
print(p)
dev.off()

# gene sets scorcards 
scorecard1 <- c("MEIS2", "GATA1", "GATA2", "RUNX1", "GYPA", "GPI", "LMO2", 
                "SPI1", "SPN", "MYC", "MYB")
scorecard1 <- rev(factor(scorecard1))

scorecard2 <- c("MEIS2", "GATA1", "GATA2", "RUNX1", "ETV6", "ITGA2B", "GYPA", 
                "HBG1", "KLF1", "LMO2", "SPI1", "MYC", "MYB")
scorecard2 <- rev(factor(scorecard2))



#### sample ordering scorecard1
# Erythroid expression
erys <- subset(res_renamed, subset = cell_type %in% 
                 c("HPC - Ery 1", "HPC - Ery 2", "HPC - Ery 3"))
erys@meta.data$sample <- factor(erys@meta.data$sample,
                                 levels = c("Euploid/wtGATA1/iPSC/D7", "Euploid/GATA1s/iPSC/D7",
                                            "T21/wtGATA1/iPSC/D7" ,"T21/GATA1s/iPSC/D7",
                                            "Euploid/wtGATA1/FL/W15", "Euploid/wtGATA1/FL/W15/D3",
                                            "Euploid/wtGATA1/FL/W22", "T21/wtGATA1/FL/W15", "T21/FL/W16", 
                                            "T21/wtGATA1/FL/W22", "T21/GATA1s/TMD/22%", "T21/GATA1s/TMD/68%"))
DotPlot(erys, features=scorecard1, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/pub/timecourse-d7-and-primary-scorecard1-ery.pdf", height = 6)

# Megakaryocyte expression
mks <- subset(res_renamed, subset = cell_type %in% 
                c("HPC - Mk 1", "HPC - Mk 2", "HPC - Mk 3"))
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
ggsave("plots/pub/timecourse-d7-and-primary-scorecard1-mk.pdf", height = 6)

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
ggsave("plots/pub/timecourse-d7-and-primary-scorecard2-ery.pdf", height = 7)

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
ggsave("plots/pub/timecourse-d7-and-primary-scorecard2-mk.pdf", height = 7)


####### combine scorecards 
scorecard <- unique(c(scorecard1, scorecard2))

#### sample ordering1
# Erythroid expression
erys <- subset(res_renamed, subset = cell_type %in% 
                 c("HPC - Ery 1", "HPC - Ery 2", "HPC - Ery 3"))
erys@meta.data$sample <- factor(erys@meta.data$sample,
                                levels = c("Euploid/wtGATA1/iPSC/D7", "Euploid/GATA1s/iPSC/D7",
                                           "T21/wtGATA1/iPSC/D7" ,"T21/GATA1s/iPSC/D7",
                                           "Euploid/wtGATA1/FL/W15", "Euploid/wtGATA1/FL/W15/D3",
                                           "Euploid/wtGATA1/FL/W22", "T21/wtGATA1/FL/W15", "T21/FL/W16", 
                                           "T21/wtGATA1/FL/W22", "T21/GATA1s/TMD/22%", "T21/GATA1s/TMD/68%"))
DotPlot(erys, features=scorecard, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/pub/timecourse-d7-and-primary-ordering1-ery.pdf", height = 7)

# Megakaryocyte expression
mks <- subset(res_renamed, subset = cell_type %in% 
                c("HPC - Mk 1", "HPC - Mk 2", "HPC - Mk 3"))

mks@meta.data$sample <- factor(mks@meta.data$sample,
                               levels = c("Euploid/wtGATA1/iPSC/D7", "Euploid/GATA1s/iPSC/D7",
                                          "T21/wtGATA1/iPSC/D7" ,"T21/GATA1s/iPSC/D7",
                                          "Euploid/wtGATA1/FL/W15", "Euploid/wtGATA1/FL/W15/D3",
                                          "Euploid/wtGATA1/FL/W22", "T21/wtGATA1/FL/W15", "T21/FL/W16", 
                                          "T21/wtGATA1/FL/W22", "T21/GATA1s/TMD/22%", "T21/GATA1s/TMD/68%"))
DotPlot(mks , features=scorecard, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/pub/timecourse-d7-and-primary-ordering1-mk.pdf", height = 7)


# Myeloid expression
myeloids <- subset(res_renamed, subset = cell_type %in% 
                c("HPC - Myeloid 1", "HPC - Myeloid 2", "HPC - Myeloid 3"))

myeloids@meta.data$sample <- factor(myeloids@meta.data$sample,
                               levels = c("Euploid/wtGATA1/iPSC/D7", "Euploid/GATA1s/iPSC/D7",
                                          "T21/wtGATA1/iPSC/D7" ,"T21/GATA1s/iPSC/D7",
                                          "Euploid/wtGATA1/FL/W15", "Euploid/wtGATA1/FL/W15/D3",
                                          "Euploid/wtGATA1/FL/W22", "T21/wtGATA1/FL/W15", "T21/FL/W16", 
                                          "T21/wtGATA1/FL/W22", "T21/GATA1s/TMD/22%", "T21/GATA1s/TMD/68%"))
DotPlot(myeloids , features=scorecard, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/pub/timecourse-d7-and-primary-ordering1-myeloid.pdf", height = 7)



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
ggsave("plots/pub/timecourse-d7-and-primary-ordering2-ery.pdf", height = 7)


# create violin plots
pdf(file = "plots/pub/timecourse-d7-and-primary-ordering2-ery-violin.pdf")
gene_list <- c("HBG1", "KLF1", "GATA1", "GATA2", "GYPA", "GPI", "LMO2", "MYC")
VlnPlot(erys, gene_list, stack = TRUE, flip = TRUE, group.by = "sample") +
  theme(legend.position = "none")
dev.off()

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
ggsave("plots/pub/timecourse-d7-and-primary-ordering2-mk.pdf", height = 7)


# create violin plots
pdf(file = "plots/pub/timecourse-d7-and-primary-ordering2-mk-violin.pdf")
gene_list <- c("ITGA2B", "GATA1", "GATA2", "RUNX1", "GPI", "LMO2", "MYC")
VlnPlot(erys, gene_list, stack = TRUE, flip = TRUE, group.by = "sample") +
  theme(legend.position = "none")
dev.off()


# Myeloid expression
myeloids@meta.data$sample <- factor(myeloids@meta.data$sample,
                               levels = c("Euploid/wtGATA1/iPSC/D7","Euploid/wtGATA1/FL/W15",
                                          "Euploid/wtGATA1/FL/W15/D3", "Euploid/wtGATA1/FL/W22", 
                                          "Euploid/GATA1s/iPSC/D7", "T21/wtGATA1/iPSC/D7", 
                                          "T21/wtGATA1/FL/W15", "T21/FL/W16", "T21/wtGATA1/FL/W22",
                                          "T21/GATA1s/iPSC/D7", "T21/GATA1s/TMD/22%", "T21/GATA1s/TMD/68%"))
DotPlot(myeloids, features=scorecard, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/pub/timecourse-d7-and-primary-ordering2-myeloid.pdf", height = 7)

# create violin plots
pdf(file = "plots/pub/timecourse-d7-and-primary-ordering2-myeloid-violin.pdf")
gene_list <- c("GATA1", "GATA2", "GPI", "LMO2", "SPI1", "SPN", "MYC", "MYB")
VlnPlot(erys, gene_list, stack = TRUE, flip = TRUE, group.by = "sample") +
  theme(legend.position = "none")
dev.off()

# expression dotPlot of key genes
hpc <- c("ITGA4", "RUNX1", "ANGPT1", "FLI1", "NFE2", "GATA2")
meg <- c("PF4", "GP9", "DIAPH1", "MMRN1", "FLNA", "DAAM1", "VCL", "WDR1",  
         "ITGB3", "ITGA2B")
ery <- c("HBZ", "KLF1", "HBE1", "GYPA", "PRDX2", "HBG2", "HBA1", "HBA2", 
         "HBG1" , "ALAS2")
myeloid <- c("CCL3", "CCL4", "MRC1", "NRP1", "SAMHD1", "CD74", "CST3")

key_genes <- c(hpc, meg, ery, myeloid)

P1 <- DotPlot(res_renamed, features = key_genes, cols = c("#f7e69c", "red"),
              dot.scale = 4)

d7Primary_data <- P1$data
d7Primary_data$condition<-"d7Primary"

plot_data <- d7Primary_data
plot_data$condition <- factor(plot_data$condition, levels = "d7Primary")
plot_data$id <- factor(plot_data$id, 
                       levels = c("HPC - unbiased", 
                                  "HPC - Mk 1", "HPC - Mk 2", "HPC - Mk 3", 
                                  "HPC - Ery 1", "HPC - Ery 2", "HPC - Ery 3", 
                                  "HPC - Myeloid 1", "HPC - Myeloid 2", "HPC - Myeloid 3"))
plot_data$features.plot <- factor(plot_data$features.plot, 
                                  levels = rev(levels(plot_data$features.plot)))
cols <- rep(c("turquoise", "palevioletred1", "mediumpurple1", "snow3"), c(7, 10, 10, 6)) 
pdf(file = "plots/pub/timecourse-d7-and-primary-umap-marker-dotplots.pdf", width = 5, height = 9)
P2 <- ggplot(plot_data, aes(x=features.plot, y=id, colour=avg.exp.scaled)) +
  geom_point(aes(size=pct.exp)) +
  scale_colour_gradient(low = "#f7e69c", high = "red") + theme_classic() + coord_flip() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.y = element_text(colour = cols), axis.title.y = element_blank())
print(P2)
dev.off()

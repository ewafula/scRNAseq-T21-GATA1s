# Annotate integrated samples and DE analysis

# Eric Wafula
# 10/20/2024


library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(scCustomize)

# load Seurat RDS Object
d7911 <- readRDS("seurat_objects/timecourse-d7911.RDS")

# get clustering results
res <- FindClusters(d7911, resolution = 0.4)

pdf(file = "plots/timecourse/timecourse-d7911-umap.pdf")
Seurat::DimPlot(res, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
dev.off()

#####  eythroid genes #######
# https://maayanlab.cloud/Harmonizome/gene_set/erythrocyte/TISSUES+Curated+Tissue+Protein+Expression+Evidence+Scores
ery_genes <- c("HBA", "HBB", "HBZ", "KLF1", "HBE1", "GYPA", "GYPB","BBS1", "DMTN", "GSR", 
               "LANCL1", "MPP1", "ATP2B4", "ATP2B1", "CYB5A","HDHD1","PLSCR1",
               "ALDH1A1","CAT","CD47", "BLVRB", "MT1E", "CD44","TSTA3", "ACHE",
               "EPB41","EPB42", "GLO1","BMI1", "GATA1", "ADD2", "PCMT1","CYB5R3",
               "ADD1", "PRDX6", "PRDX2", "PGK1", "STOM", "CD58","CTSE","SKP1",
               "HBG2", "HBA1", "HBA2", "HBG1", "SLC25A37", "ANK1", "TFRC", "CD36",
               "SPTA1", "BPGM", "ALAS2")

# add scores
res_scores <- Seurat::AddModuleScore(res, features = list(ery_genes),
                                     name="ery_enriched")
# visualize, select key genes, identify clusters
Seurat::DotPlot(res, features = ery_genes, cols = c("lightgrey", "red"))
pdf(file = "plots/timecourse/timecourse-d7911-erythroids-umap.pdf")
Seurat::FeaturePlot(res_scores, features = "ery_enriched1",
                    label = FALSE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

key_ery_genes1 <- c("HBZ", "KLF1", "HBE1", "GYPA", "BLVRB", "PRDX2", "HBG2", 
                    "HBA1", "HBA2", "HBG1" , "ALAS2")
pdf(file = "plots/timecourse/timecourse-d7911-erythroids-key-genes-umap.pdf", width = 12, height = 10)
FeaturePlot(res, features = key_ery_genes1, min.cutoff = 0, max.cutoff = 7.5)
dev.off()

Seurat::DotPlot(res, features = key_ery_genes1, cols = c("lightgrey", "red"))


#####  megakaryocytes genes #######
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
# add scores
res_scores <- Seurat::AddModuleScore(res, features = list(mk_genes),
                                     name="mk_enriched")
# visualize, select key genes, identify clusters
Seurat::DotPlot(res, features = mk_genes, cols = c("lightgrey", "red"))
pdf(file = "plots/timecourse/timecourse-d7911-megakaryocytes-umap.pdf")
Seurat::FeaturePlot(res_scores, features = "mk_enriched1", label = FALSE,
                    repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

key_mk_genes1 <- c("PF4", "GP9", "DIAPH1", "MMRN1", "FLNA", "DAAM1", "VCL", "WDR1",  
                   "LTBP1", "ITGB3", "ITGA2B", "MED12L", "PRKAR2B", "ARHGAP6")
pdf(file = "plots/timecourse/timecourse-d7911-megakaryocytes-key-genes-umap.pdf", width = 12, height = 10)
Seurat::FeaturePlot(res, features = key_mk_genes1, min.cutoff = 0, max.cutoff = 7.5)
dev.off()

# Seurat::DotPlot(res, features = c(key_mk_genes1), cols = c("lightgrey", "red"))
Seurat::DotPlot(res, features = c(key_ery_genes1, key_mk_genes1),
                cols = c("lightgrey", "red"))


#####  myeloid genes #######
# from Chris except the EMP3 and VASP from MsigDB
# https://www.gsea-msigdb.org/gsea/msigdb/cards/MA_MYELOID_DIFFERENTIATION_UP
myeloid_genes <- c("ACSL1", "AZU1", "CEBPA", "CSTA", "CXCL8", "CYBB", "CCL3",
                   "CCL4", "CCL4l2", "CLC", "MPO", "MRC1", "IRAK3", "NRP1", "IKZF2",
                   "PDE4D", "MS4A3", "RBM47", "PRTN3", "PRG2", "PRG3", "SAMHD1",
                   "S100A8", "S100A9", "VIM", "PU1", "SP1", "EMP3", "VASP", "CD74", 
                   "SRGN", "CST3", "AOAH", "SLC9A9", "CTSB", "LGMN", "FTL")

# add scores
res_scores <- Seurat::AddModuleScore(res, features = list(myeloid_genes),
                                     name="myeloid_enriched")
# visualize, select key genes, identify clusters
Seurat::DotPlot(res, features = myeloid_genes, cols = c("lightgrey", "red"))
pdf(file = "plots/timecourse/timecourse-d7911-myeloids-umap.pdf")
Seurat::FeaturePlot(res_scores, features = "myeloid_enriched1",
                    label = FALSE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

key_myeloid_genes1 <- c("CCL3", "CCL4", "MRC1", "NRP1", "RBM47", "SAMHD1", "CD74", 
                        "CST3")
pdf(file = "plots/timecourse/timecourse-d7911-myeloids-key-genes-umap.pdf", width = 10, height = 10)
FeaturePlot(res, features = key_myeloid_genes1, min.cutoff = 0, max.cutoff = 7.5)
dev.off()

# Seurat::DotPlot(res, features = c( key_myeloid_genes1), cols = c("lightgrey", "red"))
Seurat::DotPlot(res, features = c(key_ery_genes1, key_mk_genes1, key_myeloid_genes1),
                cols = c("lightgrey", "red"))
pdf(file = "plots/timecourse/timecourse-d7911-key-genes-doplot.pdf", width = 30, height = 7)
Seurat::DotPlot(res, features = c(key_ery_genes1, key_mk_genes1, key_myeloid_genes1),
                cols = c("lightgrey", "red"))
dev.off()

# Rename identity classes
res_renamed <- Seurat::RenameIdents(res, "0" = "HPC - MK bias 1", 
                                    "1" = "HPC - MK bias 2", "2" = "Ery", 
                                    "3" = "HPC - MK bias 1", "4" = "Myeloid", 
                                    "5" = "HPC - Ery bias", "6" = "Myeloid", 
                                    "7" = "MK", "8" = "MK", "9" = "HPCs",  
                                    "10" = "HPC - MK bias 2", "11" = "HPCs", 
                                    "12" = "HPCs")

Idents(res_renamed) <-  factor(Idents(res_renamed), 
                               levels = c("HPCs", "HPC - MK bias 1", "HPC - MK bias 2",
                                          "MK", "HPC - Ery bias", "Ery", "Myeloid"))

pdf(file = "plots/timecourse/timecourse-d7911-umap-renamed.pdf", width = 10)
Seurat::DimPlot(res_renamed, reduction = "umap", raster = FALSE,
                cols = c("HPCs" = "snow3", "HPC - MK bias 2" = "thistle",
                         "HPC - MK bias 1" = "mediumpurple1", "MK" = "purple4",
                         "HPC - Ery bias" = "hotpink3", "Ery" = "red4",
                         "Myeloid" = "turquoise"))
dev.off()
pdf(file = "plots/timecourse/timecourse-d7911-umap-renamed.v2.pdf")
Seurat::DimPlot(res_renamed, reduction = "umap", raster = FALSE,
                cols = c("HPCs" = "snow3", "HPC - MK bias 2" = "thistle",
                         "HPC - MK bias 1" = "mediumpurple1", "MK" = "purple4",
                         "HPC - Ery bias" = "hotpink3", "Ery" = "red4",
                         "Myeloid" = "turquoise")) + NoLegend()
dev.off()


# add cell type annotations as column in meta data
res_renamed[["cell_type"]] <- Seurat::Idents(res_renamed)

# assign Cell-Cycle Scores
res_renamed <- Seurat::CellCycleScoring(res_renamed, s.features = cc.genes$s.genes, 
                                        g2m.features = cc.genes$g2m.genes)

# save renamed seurat object
saveRDS(res_renamed, file = file.path("seurat_objects", "timecourse-d7911-renamed.RDS"))

# write cell counts to file
table(res_renamed@meta.data$Sample, res_renamed@meta.data$cell_type) %>%
  as.data.frame.matrix() %>%
  tibble::rownames_to_column(var = "Sample") %>%
  readr::write_tsv("results/timecourse/timecourse-d7911-umap-renamed-cell-counts.tsv.gz")

pdf(file = "plots/timecourse/timecourse-d7911-umap-renamed-split.pdf", width = 70)
Seurat::DimPlot(res_renamed, reduction = "umap", split.by = "Sample",
                cols = c("HPCs" = "snow3", "HPC - MK bias 2" = "thistle",
                         "HPC - MK bias 1" = "mediumpurple1", "MK" = "purple4",
                         "HPC - Ery bias" = "hotpink3", "Ery" = "red4",
                         "Myeloid" = "turquoise"))
dev.off()

# Sample cluster proportion 
cols = c("HPCs" = "snow3", "HPC - MK bias 1" = "thistle", 
         "HPC - MK bias 2" = "mediumpurple1", "MK" = "purple4",
         "HPC - Ery bias" = "hotpink3", "Ery" = "red4", "Myeloid" = "turquoise") 
tt<-prop.table(table(res_renamed$cell_type, res_renamed$Sample), margin = 2)
df<-melt(tt)
colnames(df)<-c("CellType", "Sample", "Proportion")

## D7 sample proportions
df_d7 <- df %>%  dplyr::filter(Sample %in% c("T21wtGATA1D7", "T21GATA1sD7", 
                                             "EuploidwtGATA1D7", "EuploidGATA1sD7"))
df_d7 <- df_d7 %>% dplyr::mutate(Sample = case_when(Sample == "T21wtGATA1D7" ~ "T21/wtGATA1",
                                              Sample == "T21GATA1sD7" ~ "T21/GATA1s",
                                              Sample == "EuploidwtGATA1D7" ~ "Euploid/wtGATA1",
                                              Sample == "EuploidGATA1sD7" ~ "Euploid/GATA1s"))
df_d7$CellType<-factor(df_d7$CellType, levels = rev(c("HPCs", "HPC - MK bias 1", "HPC - MK bias 2",
                                                      "MK", "HPC - Ery bias", "Ery", "Myeloid")))
df_d7$Sample<-factor(df_d7$Sample, levels = rev(c("T21/wtGATA1", "T21/GATA1s", 
                                                  "Euploid/wtGATA1", "Euploid/GATA1s")))
pdf(file = "plots/timecourse/timecourse-d7911-umap-renamed-proportions-d7.pdf")
p<-ggplot(df_d7, aes(Proportion, Sample, fill = CellType)) + geom_bar(stat = "identity") +
  scale_fill_manual(values = cols, name = "Cell Type", guide = guide_legend(reverse = TRUE)) +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.line.x = element_blank(),
        axis.line.y = element_blank(), axis.title.y = element_blank(),
        aspect.ratio = .4)
print(p)
dev.off()

## D9 sample proportions
df_d9 <- df %>%  dplyr::filter(Sample %in% c("T21wtGATA1D9", "T21GATA1sD9", 
                                             "EuploidwtGATA1D9", "EuploidGATA1sD9"))
df_d9 <- df_d9 %>% dplyr::mutate(Sample = case_when(Sample == "T21wtGATA1D9" ~ "T21/wtGATA1",
                                                    Sample == "T21GATA1sD9" ~ "T21/GATA1s",
                                                    Sample == "EuploidwtGATA1D9" ~ "Euploid/wtGATA1",
                                                    Sample == "EuploidGATA1sD9" ~ "Euploid/GATA1s"))
df_d9$CellType<-factor(df_d9$CellType, levels = rev(c("HPCs", "HPC - MK bias 1", "HPC - MK bias 2",
                                                      "MK", "HPC - Ery bias", "Ery", "Myeloid")))
df_d9$Sample<-factor(df_d9$Sample, levels = rev(c("T21/wtGATA1", "T21/GATA1s", 
                                                  "Euploid/wtGATA1", "Euploid/GATA1s")))
pdf(file = "plots/timecourse/timecourse-d7911-umap-renamed-proportions-d9.pdf")
p<-ggplot(df_d9, aes(Proportion, Sample, fill = CellType)) + geom_bar(stat = "identity") +
  scale_fill_manual(values = cols, name = "Cell Type", guide = guide_legend(reverse = TRUE)) +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.line.x = element_blank(),
        axis.line.y = element_blank(), axis.title.y = element_blank(),
        aspect.ratio = .4)
print(p)
dev.off()

## D11 sample proportions
df_d11 <- df %>%  dplyr::filter(Sample %in% c("T21wtGATA1D11", "T21GATA1sD11", 
                                             "EuploidwtGATA1D11", "EuploidGATA1sD11"))
df_d11 <- df_d11 %>% dplyr::mutate(Sample = case_when(Sample == "T21wtGATA1D11" ~ "T21/wtGATA1",
                                                    Sample == "T21GATA1sD11" ~ "T21/GATA1s",
                                                    Sample == "EuploidwtGATA1D11" ~ "Euploid/wtGATA1",
                                                    Sample == "EuploidGATA1sD11" ~ "Euploid/GATA1s"))
df_d11$CellType<-factor(df_d11$CellType, levels = rev(c("HPCs", "HPC - MK bias 1", "HPC - MK bias 2",
                                                      "MK", "HPC - Ery bias", "Ery", "Myeloid")))
df_d11$Sample<-factor(df_d11$Sample, levels = rev(c("T21/wtGATA1", "T21/GATA1s", 
                                                  "Euploid/wtGATA1", "Euploid/GATA1s")))
pdf(file = "plots/timecourse/timecourse-d7911-umap-renamed-proportions-d11.pdf")
p<-ggplot(df_d11, aes(Proportion, Sample, fill = CellType)) + geom_bar(stat = "identity") +
  scale_fill_manual(values = cols, name = "Cell Type", guide = guide_legend(reverse = TRUE)) +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.line.x = element_blank(),
        axis.line.y = element_blank(), axis.title.y = element_blank(),
        aspect.ratio = .4)
print(p)
dev.off()


# expression dotPlot of key genes
hpc <- c("ITGA4", "RUNX1", "ANGPT1", "FLI1", "NFE2", "GATA2")
meg <- c("PF4", "GP9", "DIAPH1", "MMRN1", "FLNA", "DAAM1", "VCL", "WDR1",  
         "ITGB3", "ITGA2B")
ery <- c("HBZ", "KLF1", "HBE1", "GYPA", "PRDX2", "HBG2", "HBA1", "HBA2", 
         "HBG1" , "ALAS2")
myeloid <- c("CCL3", "CCL4", "MRC1", "NRP1", "SAMHD1", "CD74", "CST3")

key_genes <- c(hpc, meg, ery, myeloid)

P1 <- DotPlot(res_renamed, features = key_genes, cols = c("#f7e69c", "red"), dot.scale = 4)

d7911_data <- P1$data
d7911_data$condition<-"d7911"

plot_data <- d7911_data
plot_data$condition <- factor(plot_data$condition, levels = "d7911")
plot_data$id <- factor(plot_data$id, 
                       levels = c("HPCs", "HPC - MK bias 2", "HPC - MK bias 1", 
                                  "MK", "HPC - Ery bias", "Ery", "Myeloid"))
plot_data$features.plot <- factor(plot_data$features.plot, 
                                        levels = rev(levels(plot_data$features.plot)))
cols <- rep(c("turquoise",  "red4", "purple4", "snow3"), c(7, 10, 10, 6)) 
pdf(file = "plots/timecourse/timecourse-d7911-umap-renamed-marker-dotplots.pdf", width = 5, height = 9)
P2 <- ggplot(plot_data, aes(x=features.plot, y=id, colour=avg.exp.scaled)) +
  geom_point(aes(size=pct.exp)) +
  scale_colour_gradient(low = "#f7e69c", high = "red") + theme_classic() + coord_flip() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = cols), axis.title.y = element_blank())
print(P2)
dev.off()

# set default assay to SCT
Seurat::DefaultAssay(res_renamed) <- "SCT"
 
####### combine scorecards 
# gene sets scorcards 
scorecard1 <- c("MEIS2", "GATA1", "GATA2", "RUNX1", "GYPA", "GPI", "LMO2", 
                "SPI1", "SPN", "MYC", "MYB")
scorecard1 <- rev(factor(scorecard1))

scorecard2 <- c("MEIS2", "GATA1", "GATA2", "RUNX1", "ETV6", "ITGA2B", "GYPA", 
                "HBG1", "KLF1", "LMO2", "SPI1", "MYC", "MYB")
scorecard2 <- rev(factor(scorecard2))
scorecard <- unique(c(scorecard1, scorecard2))

#### sample ordering1
# Erythroid expression
erys <- subset(res_renamed, subset = cell_type %in% 
                 c("HPC - Ery bias", "Ery "))
erys@meta.data$Sample <- factor(erys@meta.data$Sample,
                                levels = c("EuploidwtGATA1D7", "EuploidGATA1sD7",
                                           "T21wtGATA1D7", "T21GATA1sD7",
                                           "EuploidwtGATA1D9", "EuploidGATA1sD9",
                                           "T21wtGATA1D9", "T21GATA1sD9", 
                                           "EuploidwtGATA1D11", "EuploidGATA1sD11",
                                           "T21wtGATA1D11", "T21GATA1sD11"))
DotPlot(erys, features=scorecard, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="Sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/timecourse/timecourse-d7911-ordering1-ery.pdf", height = 7)

# Megakaryocyte expression
mks <- subset(res_renamed, subset = cell_type %in% 
                c("HPC - MK bias 1", "HPC - MK bias 2", "MK"))

mks@meta.data$Sample <- factor(mks@meta.data$Sample,
                               levels = c("EuploidwtGATA1D7", "EuploidGATA1sD7",
                                                   "T21wtGATA1D7", "T21GATA1sD7",
                                                   "EuploidwtGATA1D9", "EuploidGATA1sD9",
                                                   "T21wtGATA1D9", "T21GATA1sD9", 
                                                   "EuploidwtGATA1D11", "EuploidGATA1sD11",
                                                   "T21wtGATA1D11", "T21GATA1sD11"))
DotPlot(mks , features=scorecard, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="Sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/timecourse/timecourse-d7911-ordering1-mk.pdf", height = 7)


# Myeloid expression
myeloids <- subset(res_renamed, subset = cell_type == "Myeloid")

myeloids@meta.data$Sample <- factor(myeloids@meta.data$Sample,
                                    levels = c("EuploidwtGATA1D7", "EuploidGATA1sD7",
                                               "T21wtGATA1D7", "T21GATA1sD7",
                                               "EuploidwtGATA1D9", "EuploidGATA1sD9",
                                               "T21wtGATA1D9", "T21GATA1sD9", 
                                               "EuploidwtGATA1D11", "EuploidGATA1sD11",
                                               "T21wtGATA1D11", "T21GATA1sD11"))
DotPlot(myeloids , features=scorecard, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="Sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/timecourse/timecourse-d7911-ordering1-myeloid.pdf", height = 7)



#### sample ordering2
# Erythroid expression
erys@meta.data$Sample <- factor(erys@meta.data$Sample,
                                levels = c("EuploidwtGATA1D7", "EuploidGATA1sD7",
                                           "EuploidwtGATA1D9", "EuploidGATA1sD9",
                                           "EuploidwtGATA1D11", "EuploidGATA1sD11",
                                           "T21wtGATA1D7", "T21GATA1sD7",
                                           "T21wtGATA1D9", "T21GATA1sD9", 
                                           "T21wtGATA1D11", "T21GATA1sD11"))
DotPlot(erys, features=scorecard, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="Sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/timecourse/timecourse-d7911-ordering2-ery.pdf", height = 7)

# Megakaryocyte expression
mks@meta.data$Sample <- factor(mks@meta.data$Sample,
                               levels = c("EuploidwtGATA1D7", "EuploidGATA1sD7",
                                          "EuploidwtGATA1D9", "EuploidGATA1sD9",
                                          "EuploidwtGATA1D11", "EuploidGATA1sD11",
                                          "T21wtGATA1D7", "T21GATA1sD7",
                                          "T21wtGATA1D9", "T21GATA1sD9", 
                                          "T21wtGATA1D11", "T21GATA1sD11"))
DotPlot(mks, features=scorecard, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="Sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/timecourse/timecourse-d7911-ordering2-mk.pdf", height = 7)


# Myeloid expression
myeloids@meta.data$Sample <- factor(myeloids@meta.data$Sample,
                                    levels = c("EuploidwtGATA1D7", "EuploidGATA1sD7",
                                               "EuploidwtGATA1D9", "EuploidGATA1sD9",
                                               "EuploidwtGATA1D11", "EuploidGATA1sD11",
                                               "T21wtGATA1D7", "T21GATA1sD7",
                                               "T21wtGATA1D9", "T21GATA1sD9", 
                                               "T21wtGATA1D11", "T21GATA1sD11"))
DotPlot(myeloids, features=scorecard, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="Sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/timecourse/timecourse-d7911-ordering2-myeloid.pdf", height = 7)


#### sample ordering3
# Erythroid expression
erys@meta.data$Sample <- factor(erys@meta.data$Sample,
                                levels = c("EuploidwtGATA1D7", "EuploidwtGATA1D9", "EuploidwtGATA1D11",
                                           "EuploidGATA1sD7", "EuploidGATA1sD9", "EuploidGATA1sD11",
                                           "T21wtGATA1D7", "T21wtGATA1D9", "T21wtGATA1D11", 
                                           "T21GATA1sD7", "T21GATA1sD9", "T21GATA1sD11"))

DotPlot(erys, features=scorecard, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="Sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/timecourse/timecourse-d7911-ordering3-ery.pdf", height = 7)

# create violin plots
pdf(file = "plots/pub/timecourse-d7911-ordering3-ery-violin.pdf")
gene_list <- c("HBG1", "KLF1", "MEIS2", "GATA1", "GATA2", "RUNX1", "GYPA")
VlnPlot(erys, gene_list, stack = TRUE, flip = TRUE, group.by = "Sample") +
  theme(legend.position = "none")
dev.off()

# Megakaryocyte expression
mks@meta.data$Sample <- factor(mks@meta.data$Sample,
                               levels = c("EuploidwtGATA1D7", "EuploidwtGATA1D9", "EuploidwtGATA1D11",
                                          "EuploidGATA1sD7", "EuploidGATA1sD9", "EuploidGATA1sD11",
                                          "T21wtGATA1D7", "T21wtGATA1D9", "T21wtGATA1D11", 
                                          "T21GATA1sD7", "T21GATA1sD9", "T21GATA1sD11"))
DotPlot(mks, features=scorecard, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="Sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/timecourse/timecourse-d7911-ordering3-mk.pdf", height = 7)

# create violin plots
pdf(file = "plots/pub/timecourse-d7911-ordering3-mk-violin.pdf")
gene_list <- c("ETV6", "ITGA2B", "MEIS2", "GATA1", "GATA2", "RUNX1")
VlnPlot(mks, gene_list, stack = TRUE, flip = TRUE, group.by = "Sample") +
  theme(legend.position = "none")
dev.off()

# Myeloid expression
myeloids@meta.data$Sample <- factor(myeloids@meta.data$Sample,
                                    levels = c("EuploidwtGATA1D7", "EuploidwtGATA1D9", "EuploidwtGATA1D11",
                                                "EuploidGATA1sD7", "EuploidGATA1sD9", "EuploidGATA1sD11",
                                                "T21wtGATA1D7", "T21wtGATA1D9", "T21wtGATA1D11", 
                                                "T21GATA1sD7", "T21GATA1sD9", "T21GATA1sD11"))
DotPlot(myeloids, features=scorecard, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="Sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/timecourse/timecourse-d7911-ordering3-myeloid.pdf", height = 7)

# create violin plots
pdf(file = "plots/pub/timecourse-d7911-ordering3-myeloid-violin.pdf")
gene_list <- c("MEIS2", "GATA1", "GATA2", "RUNX1", "LMO2", "SPI1")
VlnPlot(myeloids, gene_list, stack = TRUE, flip = TRUE, group.by = "Sample") +
  theme(legend.position = "none")
dev.off()

# HPC expression
hpcs <- subset(res_renamed, subset = cell_type == "HPCs")
hpcs@meta.data$Sample <- factor(hpcs@meta.data$Sample,
                                    levels = c("EuploidwtGATA1D7", "EuploidwtGATA1D9", "EuploidwtGATA1D11",
                                               "EuploidGATA1sD7", "EuploidGATA1sD9", "EuploidGATA1sD11",
                                               "T21wtGATA1D7", "T21wtGATA1D9", "T21wtGATA1D11", 
                                               "T21GATA1sD7", "T21GATA1sD9", "T21GATA1sD11"))
DotPlot(myeloids, features=scorecard, cols=c("#f7e69c","red3"), dot.scale = 4,
        group.by="Sample") + coord_flip() + 
  theme(axis.text.x=element_text(angle=55, hjust=1))
ggsave("plots/timecourse/timecourse-d7911-ordering3-hpcs.pdf", height = 7)

# create violin plots
pdf(file = "plots/pub/timecourse-d7911-ordering3-hpc-violin.pdf")
gene_list <- c("ITGA2B", "KLF1", "MEIS2", "GATA1", "GATA2", "RUNX1", "LMO2")
VlnPlot(hpcs, gene_list, stack = TRUE, flip = TRUE, group.by = "Sample") +
  theme(legend.position = "none")
dev.off()

# reset default assay to integrated
Seurat::DefaultAssay(res_renamed) <- "integrated"

# update sample names to plot sample-specific cell-type UMAPs 
res_renamed_tmp <- res_renamed
res_renamed_tmp@meta.data <- res_renamed_tmp@meta.data %>% 
  dplyr::mutate(Sample = 
                  case_when(Sample == "T21wtGATA1D7" ~ "T21/wtGATA1/D7", 
                            Sample == "T21GATA1sD7" ~ "T21/GATA1s/D7", 
                            Sample == "EuploidwtGATA1D7" ~ "Euploid/wtGATA1/D7", 
                            Sample == "EuploidGATA1sD7" ~ "Euploid/GATA1s/D7" , 
                            Sample == "T21wtGATA1D9" ~ "T21/wtGATA1/D9",
                            Sample == "T21GATA1sD9" ~ "T21/GATA1s/D9",
                            Sample == "EuploidwtGATA1D9" ~ "Euploid/wtGATA1/D9",
                            Sample == "EuploidGATA1sD9" ~ "Euploid/GATA1s/D9",
                            Sample == "T21wtGATA1D11"  ~ "T21/wtGATA1/D11" ,
                            Sample == "T21GATA1sD11"  ~ "T21/GATA1s/D11" ,
                            Sample == "EuploidwtGATA1D11" ~ "Euploid/wtGATA1/D11",
                            Sample == "EuploidGATA1sD11" ~ "Euploid/GATA1s/D11"))

# recorder updated sample names
res_renamed_tmp@meta.data$Sample <- 
  factor(res_renamed_tmp@meta.data$Sample, 
         levels = c("Euploid/wtGATA1/D7", "Euploid/wtGATA1/D9", "Euploid/wtGATA1/D11",
                    "Euploid/GATA1s/D7", "Euploid/GATA1s/D9", "Euploid/GATA1s/D11",
                    "T21/wtGATA1/D7", "T21/wtGATA1/D9", "T21/wtGATA1/D11",
                    "T21/GATA1s/D7", "T21/GATA1s/D9", "T21/GATA1s/D11"))

# plot sample-specific cell-type UMAPs
cols = c("HPCs" = "snow3", "HPC - MK bias 1" = "thistle", 
         "HPC - MK bias 2" = "mediumpurple1", "MK" = "purple4",
         "HPC - Ery bias" = "hotpink3", "Ery" = "red4", "Myeloid" = "turquoise")
scCustomize::DimPlot_scCustom(res_renamed_tmp, reduction = "umap", group.by = "cell_type", 
                              label = FALSE,  colors_use = cols, split.by = "Sample",
                              num_columns = 3, raster = FALSE, split_seurat = TRUE,
                              figure_plot = TRUE) 
ggsave("plots/timecourse/timecourse-d7911-umap-renamed-split-v2.pdf",
       width = 38, height = 48)


##### Cell cycle phases plot
# plot MK cell cycle phase
cc_phases <- res_renamed@meta.data %>% dplyr::select(Sample, cell_type, Phase) %>%
  dplyr::mutate(cell_type = case_when(cell_type == "HPCs" ~ "HPC", 
                                      cell_type == "HPC - MK bias 1" ~ "MK1", 
                                      cell_type == "HPC - MK bias 2" ~ "MK2",
                                      cell_type == "MK" ~ "MK3",
                                      cell_type == "HPC - Ery bias" ~ "Erythroid1",
                                      cell_type == "Ery" ~ "Erythroid2",
                                      TRUE ~ cell_type)) %>% 
  dplyr::mutate(genotype = case_when(Sample %in% 
                                       c("EuploidwtGATA1D7", "EuploidwtGATA1D9",
                                         "EuploidwtGATA1D11" ) ~ "Euploid/wtGATA1",
                                     Sample %in% 
                                       c("EuploidGATA1sD7", "EuploidGATA1sD9",
                                         "EuploidGATA1sD11" ) ~ "Euploid/GATA1s",
                                     Sample %in% 
                                       c("T21wtGATA1D7", "T21wtGATA1D9",
                                         "T21wtGATA1D11" ) ~ "T21/wtGATA1",
                                     Sample %in% 
                                       c("T21GATA1sD7", "T21GATA1sD9",
                                         "T21GATA1sD11" ) ~ "T21/GATA1s")) %>% 
  dplyr::select(-Sample)


# HPC
HPC <- cc_phases %>% dplyr::filter(cell_type == "HPC") %>% 
  dplyr::select(-cell_type)
HPC <- prop.table(table(HPC$Phase, HPC$genotype), 
                  margin = 2)
HPC <- HPC %>% as.data.frame() %>% 
  dplyr::rename(Phase = Var1, genotype = Var2, Proportion = Freq) %>% 
  dplyr::mutate(cell_type = "HPC")

# MK1
MK1 <- cc_phases %>% dplyr::filter(cell_type == "MK1") %>% 
  dplyr::select(-cell_type)
MK1 <- prop.table(table(MK1$Phase, MK1$genotype), 
                  margin = 2)
MK1 <- MK1 %>% as.data.frame() %>% 
  dplyr::rename(Phase = Var1, genotype = Var2, Proportion = Freq) %>% 
  dplyr::mutate(cell_type = "MK1")

# MK2
MK2 <- cc_phases %>% dplyr::filter(cell_type == "MK2") %>% 
  dplyr::select(-cell_type)
MK2 <- prop.table(table(MK2$Phase, MK2$genotype), 
                  margin = 2)
MK2 <- MK2 %>% as.data.frame() %>% 
  dplyr::rename(Phase = Var1, genotype = Var2, Proportion = Freq) %>% 
  dplyr::mutate(cell_type = "MK2")

# MK3
MK3 <- cc_phases %>% dplyr::filter(cell_type == "MK3") %>% 
  dplyr::select(-cell_type)
MK3 <- prop.table(table(MK3$Phase, MK3$genotype), 
                  margin = 2)
MK3 <- MK3 %>% as.data.frame() %>% 
  dplyr::rename(Phase = Var1, genotype = Var2, Proportion = Freq) %>% 
  dplyr::mutate(cell_type = "MK3")

# Erythroid1
Erythroid1 <- cc_phases %>% dplyr::filter(cell_type == "Erythroid1") %>% 
  dplyr::select(-cell_type)
Erythroid1 <- prop.table(table(Erythroid1$Phase, Erythroid1$genotype), 
                         margin = 2)
Erythroid1 <- Erythroid1 %>% as.data.frame() %>% 
  dplyr::rename(Phase = Var1, genotype = Var2, Proportion = Freq) %>% 
  dplyr::mutate(cell_type = "Erythroid1")

# Erythroid2
Erythroid2 <- cc_phases %>% dplyr::filter(cell_type == "Erythroid2") %>% 
  dplyr::select(-cell_type)
Erythroid2 <- prop.table(table(Erythroid2$Phase, Erythroid2$genotype), 
                         margin = 2)
Erythroid2 <- Erythroid2 %>% as.data.frame() %>% 
  dplyr::rename(Phase = Var1, genotype = Var2, Proportion = Freq) %>% 
  dplyr::mutate(cell_type = "Erythroid2")

# Myeloid
Myeloid <- cc_phases %>% dplyr::filter(cell_type == "Myeloid") %>% 
  dplyr::select(-cell_type)
Myeloid <- prop.table(table(Myeloid$Phase, Myeloid$genotype), 
                      margin = 2)
Myeloid <- Myeloid %>% as.data.frame() %>% 
  dplyr::rename(Phase = Var1, genotype = Var2, Proportion = Freq) %>% 
  dplyr::mutate(cell_type = "Myeloid")

cc_phases <-  dplyr::bind_rows(HPC, MK1, MK2, MK3, Erythroid1, Erythroid2, 
                               Myeloid)


cols = c("G2M" = "gold",  "S" = "forestgreen", "G1" = "firebrick")
cc_phases$cell_type <- factor(cc_phases$cell_type, 
                              levels = c("HPC", "MK1", "MK2", "MK3", "Erythroid1",
                                         "Erythroid2", "Myeloid"))
cc_phases$Phase <- factor(cc_phases$Phase, levels = rev(c("G2M", "S", "G1")))
cc_phases$genotype <- factor(cc_phases$genotype, 
                             levels = c("Euploid/wtGATA1", "Euploid/GATA1s",
                                        "T21/wtGATA1", "T21/GATA1s"))

pdf(file = "plots/pub/timecourse-d7911-umap-renamed-proportions-cell-cycle-phases.pdf",
    width = 12, height = 8)
ggplot2::ggplot(cc_phases, aes(genotype, Proportion, fill = Phase)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cols, name = "Phase", guide = guide_legend(reverse = TRUE)) +
  theme_classic() +
  theme(axis.text = element_text(size = 10), axis.line.x = element_blank(),
        axis.line.y = element_blank(), axis.title.x = element_blank(), 
        aspect.ratio = 1, axis.text.x = element_text(angle = 65, hjust = 1)) +
  facet_wrap(~cell_type, ncol = 4)
dev.off()


# Find Markers - DE between samples cell type with EuploidwtGATA1 as baseline reference
res_renamed_findmarkers <- res_renamed
res_renamed_findmarkers$cell_type_sample <- paste(gsub(" |/", "_", Idents(res_renamed_findmarkers)),
                                                  res_renamed_findmarkers$Sample, sep = "_")
Idents(res_renamed_findmarkers) <- "cell_type_sample"
cell_type_samples <- unique(res_renamed_findmarkers$cell_type_sample)
cell_type_samples <- cell_type_samples[!grepl("^NA_", cell_type_samples)]
res_renamed_findmarkers <- Seurat::PrepSCTFindMarkers(res_renamed_findmarkers, assay = "SCT")

# Find Markers - T21/wtGATA1 vs T21/GATA1s comparisons for D7 and D11 
# with T21/wtGATA1 as baseline reference
# D7 samples
cell_type_samples_d7 <- cell_type_samples[grepl("D7", cell_type_samples)]
cell_type_samples_d7 <- cell_type_samples_d7[grepl("T21", cell_type_samples_d7)]
cell_type_samples_d7 <- cell_type_samples_d7[!grepl("T21wtGATA1D7", cell_type_samples_d7)]
for (cell_type_sample in cell_type_samples_d7) {
  cell_type <- strsplit(cell_type_sample, split = "_")[[1]]
  control_sample <- paste(c(head(cell_type, -1), "T21wtGATA1D7"), collapse = "_")
  Seurat::FindMarkers(res_renamed_findmarkers, ident.1 = cell_type_sample, ident.2 = control_sample, 
                      min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT") %>% 
    tibble::rownames_to_column(var = "gene_symbol") %>% 
    readr::write_tsv(paste0("results/timecourse/T21GATA1s_vs_T21wtGATA1/D7/timecourse_d7911_d7_", 
                            cell_type_sample, "_dge.tsv.gz"))
}

# D11 samples
cell_type_samples_d11 <- cell_type_samples[grepl("D11", cell_type_samples)]
cell_type_samples_d11 <- cell_type_samples_d11[grepl("T21", cell_type_samples_d11)]
cell_type_samples_d11 <- cell_type_samples_d11[!grepl("T21wtGATA1D11", cell_type_samples_d11)]
for (cell_type_sample in cell_type_samples_d11) {
  cell_type <- strsplit(cell_type_sample, split = "_")[[1]]
  control_sample <- paste(c(head(cell_type, -1), "T21wtGATA1D11"), collapse = "_")
  Seurat::FindMarkers(res_renamed_findmarkers, ident.1 = cell_type_sample, ident.2 = control_sample, 
                      min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT") %>% 
    tibble::rownames_to_column(var = "gene_symbol") %>% 
    readr::write_tsv(paste0("results/timecourse/T21GATA1s_vs_T21wtGATA1/D11/timecourse_d7911_d11_", 
                            cell_type_sample, "_dge.tsv.gz"))
}


# Find Markers - Euploid/wtGATA1 vs Euploid/GATA1s comparisons for D7 and D11 
# with Euploid/wtGATA1 as baseline reference
# D7 samples
cell_type_samples_d7 <- cell_type_samples[grepl("D7", cell_type_samples)]
cell_type_samples_d7 <- cell_type_samples_d7[grepl("Euploid", cell_type_samples_d7)]
cell_type_samples_d7 <- cell_type_samples_d7[!grepl("EuploidwtGATA1D7", cell_type_samples_d7)]
for (cell_type_sample in cell_type_samples_d7) {
  cell_type <- strsplit(cell_type_sample, split = "_")[[1]]
  control_sample <- paste(c(head(cell_type, -1), "EuploidwtGATA1D7"), collapse = "_")
  Seurat::FindMarkers(res_renamed_findmarkers, ident.1 = cell_type_sample, ident.2 = control_sample, 
                      min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT") %>% 
    tibble::rownames_to_column(var = "gene_symbol") %>% 
    readr::write_tsv(paste0("results/timecourse/EuploidGATA1s_vs_EuploidwtGATA1/D7/timecourse_d7911_d7_", 
                            cell_type_sample, "_dge.tsv.gz"))
}

# D11 samples
cell_type_samples_d11 <- cell_type_samples[grepl("D11", cell_type_samples)]
cell_type_samples_d11 <- cell_type_samples_d11[grepl("Euploid", cell_type_samples_d11)]
cell_type_samples_d11 <- cell_type_samples_d11[!grepl("EuploidwtGATA1D11", cell_type_samples_d11)]
for (cell_type_sample in cell_type_samples_d11) {
  cell_type <- strsplit(cell_type_sample, split = "_")[[1]]
  control_sample <- paste(c(head(cell_type, -1), "EuploidwtGATA1D11"), collapse = "_")
  Seurat::FindMarkers(res_renamed_findmarkers, ident.1 = cell_type_sample, ident.2 = control_sample, 
                      min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT") %>% 
    tibble::rownames_to_column(var = "gene_symbol") %>% 
    readr::write_tsv(paste0("results/timecourse/EuploidGATA1s_vs_EuploidwtGATA1/D11/timecourse_d7911_d11_", 
                            cell_type_sample, "_dge.tsv.gz"))
}


# Find Markers - T21/wtGATA1 vs Euploid/wtGATA1 comparisons for D7 and D11 
# with Euploid/wtGATA1 as baseline reference
# D7 samples
cell_type_samples_d7 <- cell_type_samples[grepl("D7", cell_type_samples)]
cell_type_samples_d7 <- cell_type_samples_d7[grepl("T21", cell_type_samples_d7)]
cell_type_samples_d7 <- cell_type_samples_d7[!grepl("T21GATA1sD7", cell_type_samples_d7)]
for (cell_type_sample in cell_type_samples_d7) {
  cell_type <- strsplit(cell_type_sample, split = "_")[[1]]
  control_sample <- paste(c(head(cell_type, -1), "EuploidwtGATA1D7"), collapse = "_")
  Seurat::FindMarkers(res_renamed_findmarkers, ident.1 = cell_type_sample, ident.2 = control_sample, 
                      min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT") %>% 
    tibble::rownames_to_column(var = "gene_symbol") %>% 
    readr::write_tsv(paste0("results/timecourse/T21wtGATA1_vs_EuploidwtGATA1/D7/timecourse_d7911_d7_", 
                            cell_type_sample, "_dge.tsv.gz"))
}

# D11 samples
cell_type_samples_d11 <- cell_type_samples[grepl("D11", cell_type_samples)]
cell_type_samples_d11 <- cell_type_samples_d11[grepl("T21", cell_type_samples_d11)]
cell_type_samples_d11 <- cell_type_samples_d11[!grepl("T21GATA1sD11", cell_type_samples_d11)]
for (cell_type_sample in cell_type_samples_d11) {
  cell_type <- strsplit(cell_type_sample, split = "_")[[1]]
  control_sample <- paste(c(head(cell_type, -1), "EuploidwtGATA1D11"), collapse = "_")
  Seurat::FindMarkers(res_renamed_findmarkers, ident.1 = cell_type_sample, ident.2 = control_sample, 
                      min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT") %>% 
    tibble::rownames_to_column(var = "gene_symbol") %>% 
    readr::write_tsv(paste0("results/timecourse/T21wtGATA1_vs_EuploidwtGATA1/D11/timecourse_d7911_d11_", 
                            cell_type_sample, "_dge.tsv.gz"))
}


# GSEA analysis of integrated and annotated ssample

# Eric Wafula
# 10/20/2024



# Load libraries
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(clusterProfiler)
library(msigdbr)
library(fgsea)
library(org.Hs.eg.db)

set.seed(123)

# Retrieve MsigDB Hallmark gene sets
gene_set <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)


### T21GATAsD7 - T21wtGATA1D7 as baseline reference ####

# list dge sample cell type files
dge_files <- list.files(path = "results/timecourse/T21GATA1s_vs_T21wtGATA1/D7", 
                        pattern = "timecourse_d7911_d7_")

# run gsea on each dge file
for (dge_file in dge_files) {
  # load dge results
  dge_results <- readr::read_tsv(file.path("results/timecourse/T21GATA1s_vs_T21wtGATA1/D7", 
                                           dge_file))
  
  # pre-rank lists of fold change by sorting
  lfc_rank <- dge_results %>% 
    dplyr::pull(avg_log2FC, name = gene_symbol)
  lfc_rank <- sort(lfc_rank, decreasing = TRUE)
  
  # run gsea
  # gsea_results <- fgsea::fgsea(pathways = gene_set, stats = lfc_rank) %>%
  gsea_results <- fgsea::fgsea(pathways = gene_set, stats = lfc_rank, 
                               nPermSimple = 100000) %>% 
    tibble::as_tibble() %>% 
    dplyr::filter(padj <= 0.05) %>% 
    dplyr::arrange(dplyr::desc(NES))
  
  # write gsea results to file
  output_file <- paste0(strsplit(dge_file, 
                                split = "_dge")[[1]][1], "_gsea.tsv.gz")
  gsea_results %>%  readr::write_tsv(file.path("results/timecourse/T21GATA1s_vs_T21wtGATA1/D7", 
                                               "GSEA", output_file))
}

# list gsea sample cell type files
gsea_dir <- file.path("results/timecourse/T21GATA1s_vs_T21wtGATA1/D7", "GSEA")
gsea_files <- list.files(path = gsea_dir, pattern = "_gsea.tsv.gz")

# initialize list for gsea results
gsea_results_list <-list()

for (gsea_file in gsea_files) {
  # load dge results
  contrast <- strsplit(gsea_file, split = "_gsea")[[1]][1] 
  contrast <- gsub("timecourse_d7911_d7_", "", contrast)
  gsea_results <- readr::read_tsv(file.path(gsea_dir, gsea_file)) %>% 
    dplyr::select(pathway, padj, NES) %>% 
    dplyr::rename(GeneSet = pathway) %>% 
    dplyr::mutate(GeneSet = gsub("HALLMARK_", "", GeneSet)) %>% 
    dplyr::mutate(GeneSet = gsub("_", " ", GeneSet)) %>%
    dplyr::mutate(Group = contrast) %>% 
    dplyr::mutate(Group = case_when(Group == "HPCs_T21GATA1sD7" ~ "HPCs",
                                    Group == "HPC_-_MK_bias_2_T21GATA1sD7" ~ "HPC - MK bias 2",
                                    Group == "HPC_-_MK_bias_1_T21GATA1sD7" ~ "HPC - MK bias 1",
                                    Group == "MK_T21GATA1sD7" ~ "MK",
                                    Group == "Myeloid_T21GATA1sD7" ~ "Myeloid")) %>%
    dplyr::mutate(Enriched.group = case_when(NES > 0 ~ "1", NES < 0 ~ "2"))
  gsea_results_list <- append(gsea_results_list, list(gsea_results))
}

# merge gsea results
gsea_results_merged <- bind_rows(gsea_results_list)
gsea_results_merged$Group <- factor(gsea_results_merged$Group, 
                                    levels = c("HPCs", "HPC - MK bias 2", 
                                               "HPC - MK bias 1", "MK", "Myeloid"))

# dotplot
p <- ggplot(gsea_results_merged, aes(Group, GeneSet)) + 
  geom_point(aes(fill = Enriched.group, alpha = padj, shape =21, size = abs(NES))) + 
  geom_point(aes(color = Enriched.group, shape =21, size = abs(NES)))  +
  scale_shape_identity() 

p <- p + theme_bw() + scale_color_manual(values = c( "1" = "forestgreen", "2" = "royalblue"), 
                                         labels = c("1" = "T21/GATA1s", "2" = "T21/wtGATA1")) +
  scale_fill_manual(values = c("forestgreen", "royalblue")) + scale_alpha_continuous(range = c(1,0.1)) +
  scale_size(range = c(1,5)) + scale_x_discrete(position = "top") +
  theme(axis.text.y = element_text(hjust = 1, size = 8)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 40, hjust = 0, size = 8)) +
  # theme(axis.text.x = element_blank()) +
  # theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  labs(alpha = "Adjusted P value", size = "NES abs. value", color = "Enriched samples") +
  theme(plot.margin = margin(0.2, 2, 0.2, 0.2, "cm"))
#theme(text = element_text(face = "bold"))
pdf(file = "plots/timecourse/timecourse-d7911-dotplots-hallmark-T21GATA1s-vs-T21wtGATA1-d7.pdf")  
p <- p + guides(fill = FALSE) + theme(legend.position = "right", legend.box = "vertical")
p <- p + guides(color = guide_legend(order = 1), size  = guide_legend(order = 2), alpha = guide_legend(order = 2))
print(p)
dev.off()


### T21GATAsD11 - T21wtGATA1D11 as baseline reference ####

# list dge sample cell type files
dge_files <- list.files(path = "results/timecourse/T21GATA1s_vs_T21wtGATA1/D11", 
                        pattern = "timecourse_d7911_d11_")

# run gsea on each dge file
for (dge_file in dge_files) {
  # load dge results
  dge_results <- readr::read_tsv(file.path("results/timecourse/T21GATA1s_vs_T21wtGATA1/D11", 
                                           dge_file))
  
  # pre-rank lists of fold change by sorting
  lfc_rank <- dge_results %>% 
    dplyr::pull(avg_log2FC, name = gene_symbol)
  lfc_rank <- sort(lfc_rank, decreasing = TRUE)
  
  # run gsea
  # gsea_results <- fgsea::fgsea(pathways = gene_set, stats = lfc_rank) %>%
  gsea_results <- fgsea::fgsea(pathways = gene_set, stats = lfc_rank, 
                               nPermSimple = 100000) %>% 
    tibble::as_tibble() %>% 
    dplyr::filter(padj <= 0.05) %>% 
    dplyr::arrange(dplyr::desc(NES))
  
  # write gsea results to file
  output_file <- paste0(strsplit(dge_file, 
                                 split = "_dge")[[1]][1], "_gsea.tsv.gz")
  gsea_results %>%  readr::write_tsv(file.path("results/timecourse/T21GATA1s_vs_T21wtGATA1/D11", 
                                               "GSEA", output_file))
}

# list gsea sample cell type files
gsea_dir <- file.path("results/timecourse/T21GATA1s_vs_T21wtGATA1/D11", "GSEA")
gsea_files <- list.files(path = gsea_dir, pattern = "_gsea.tsv.gz")

# initialize list for gsea results
gsea_results_list <-list()

for (gsea_file in gsea_files) {
  # load dge results
  contrast <- strsplit(gsea_file, split = "_gsea")[[1]][1] 
  contrast <- gsub("timecourse_d7911_d11_", "", contrast)
  gsea_results <- readr::read_tsv(file.path(gsea_dir, gsea_file)) %>% 
    dplyr::select(pathway, padj, NES) %>% 
    dplyr::rename(GeneSet = pathway) %>% 
    dplyr::mutate(GeneSet = gsub("HALLMARK_", "", GeneSet)) %>% 
    dplyr::mutate(GeneSet = gsub("_", " ", GeneSet)) %>%
    dplyr::mutate(Group = contrast) %>% 
    dplyr::mutate(Group = case_when(Group == "HPC_-_MK_bias_2_T21GATA1sD11" ~ "HPC - MK bias 2",
                                    Group == "HPC_-_MK_bias_1_T21GATA1sD11" ~ "HPC - MK bias 1",
                                    Group == "MK_T21GATA1sD11" ~ "MK",
                                    Group == "HPC_-_Ery_bias_T21GATA1sD11" ~ "HPC - Ery bias",
                                    Group == "Ery_T21GATA1sD11" ~ "Ery")) %>%
    dplyr::mutate(Enriched.group = case_when(NES > 0 ~ "1", NES < 0 ~ "2"))
  gsea_results_list <- append(gsea_results_list, list(gsea_results))
}

# merge gsea results
gsea_results_merged <- bind_rows(gsea_results_list)
gsea_results_merged$Group <- factor(gsea_results_merged$Group, 
                                    levels = c("HPC - MK bias 2", "HPC - MK bias 1", 
                                               "MK", "HPC - Ery bias", "Ery"))

# dotplot
p <- ggplot(gsea_results_merged, aes(Group, GeneSet)) + 
  geom_point(aes(fill = Enriched.group, alpha = padj, shape =21, size = abs(NES))) + 
  geom_point(aes(color = Enriched.group, shape =21, size = abs(NES)))  +
  scale_shape_identity() 

p <- p + theme_bw() + scale_color_manual(values = c( "1" = "forestgreen", "2" = "royalblue"), 
                                         labels = c("1" = "T21/GATA1s", "2" = "T21/wtGATA1")) +
  scale_fill_manual(values = c("forestgreen", "royalblue")) + scale_alpha_continuous(range = c(1,0.1)) +
  scale_size(range = c(1,5)) + scale_x_discrete(position = "top") +
  theme(axis.text.y = element_text(hjust = 1, size = 8)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 40, hjust = 0, size = 8)) +
  # theme(axis.text.x = element_blank()) +
  # theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  labs(alpha = "Adjusted P value", size = "NES abs. value", color = "Enriched samples") +
  theme(plot.margin = margin(0.2, 2, 0.2, 0.2, "cm"))
#theme(text = element_text(face = "bold"))
pdf(file = "plots/timecourse/timecourse-d7911-dotplots-hallmark-T21GATA1s-vs-T21wtGATA1-d11.pdf")  
p <- p + guides(fill = FALSE) + theme(legend.position = "right", legend.box = "vertical")
p <- p + guides(color = guide_legend(order = 1), size  = guide_legend(order = 2), alpha = guide_legend(order = 2))
print(p)
dev.off()


### EuploidGATAsD7 - EuploidwtGATA1D7 as baseline reference ####

# list dge sample cell type files
dge_files <- list.files(path = "results/timecourse/EuploidGATA1s_vs_EuploidwtGATA1/D7", 
                        pattern = "timecourse_d7911_d7_")

# run gsea on each dge file
for (dge_file in dge_files) {
  # load dge results
  dge_results <- readr::read_tsv(file.path("results/timecourse/EuploidGATA1s_vs_EuploidwtGATA1/D7", 
                                           dge_file))
  
  # pre-rank lists of fold change by sorting
  lfc_rank <- dge_results %>% 
    dplyr::pull(avg_log2FC, name = gene_symbol)
  lfc_rank <- sort(lfc_rank, decreasing = TRUE)
  
  # run gsea
  # gsea_results <- fgsea::fgsea(pathways = gene_set, stats = lfc_rank) %>%
  gsea_results <- fgsea::fgsea(pathways = gene_set, stats = lfc_rank, 
                               nPermSimple = 100000) %>% 
    tibble::as_tibble() %>% 
    dplyr::filter(padj <= 0.05) %>% 
    dplyr::arrange(dplyr::desc(NES))
  
  # write gsea results to file
  output_file <- paste0(strsplit(dge_file, 
                                 split = "_dge")[[1]][1], "_gsea.tsv.gz")
  gsea_results %>%  readr::write_tsv(file.path("results/timecourse/EuploidGATA1s_vs_EuploidwtGATA1/D7", 
                                               "GSEA", output_file))
}

# list gsea sample cell type files
gsea_dir <- file.path("results/timecourse/EuploidGATA1s_vs_EuploidwtGATA1/D7", "GSEA")
gsea_files <- list.files(path = gsea_dir, pattern = "_gsea.tsv.gz")

# initialize list for gsea results
gsea_results_list <-list()

for (gsea_file in gsea_files) {
  # load dge results
  contrast <- strsplit(gsea_file, split = "_gsea")[[1]][1] 
  contrast <- gsub("timecourse_d7911_d7_", "", contrast)
  gsea_results <- readr::read_tsv(file.path(gsea_dir, gsea_file)) %>% 
    dplyr::select(pathway, padj, NES) %>% 
    dplyr::rename(GeneSet = pathway) %>% 
    dplyr::mutate(GeneSet = gsub("HALLMARK_", "", GeneSet)) %>% 
    dplyr::mutate(GeneSet = gsub("_", " ", GeneSet)) %>%
    dplyr::mutate(Group = contrast) %>% 
    dplyr::mutate(Group = case_when(Group == "HPCs_EuploidGATA1sD7" ~ "HPCs",
                                    Group == "HPC_-_MK_bias_2_EuploidGATA1sD7" ~ "HPC - MK bias 2",
                                    Group == "HPC_-_MK_bias_1_EuploidGATA1sD7" ~ "HPC - MK bias 1",
                                    Group == "MK_EuploidGATA1sD7" ~ "MK",
                                    Group == "HPC_-_Ery_bias_EuploidGATA1sD7" ~ "HPC - Ery bias",
                                    Group == "Ery_EuploidGATA1sD7" ~ "Ery",
                                    Group == "Myeloid_EuploidGATA1sD7" ~ "Myeloid")) %>%
    dplyr::mutate(Enriched.group = case_when(NES > 0 ~ "1", NES < 0 ~ "2"))
  gsea_results_list <- append(gsea_results_list, list(gsea_results))
}

# merge gsea results
gsea_results_merged <- bind_rows(gsea_results_list)
gsea_results_merged$Group <- factor(gsea_results_merged$Group, 
                                    levels = c("HPCs", "HPC - MK bias 2", "HPC - MK bias 1", 
                                               "MK", "HPC - Ery bias", "Ery", "Myeloid"))

# dotplot
p <- ggplot(gsea_results_merged, aes(Group, GeneSet)) + 
  geom_point(aes(fill = Enriched.group, alpha = padj, shape =21, size = abs(NES))) + 
  geom_point(aes(color = Enriched.group, shape =21, size = abs(NES)))  +
  scale_shape_identity() 

p <- p + theme_bw() + scale_color_manual(values = c( "1" = "forestgreen", "2" = "royalblue"), 
                                         labels = c("1" = "Euploid/GATA1s", "2" = "Euploid/wtGATA1")) +
  scale_fill_manual(values = c("forestgreen", "royalblue")) + scale_alpha_continuous(range = c(1,0.1)) +
  scale_size(range = c(1,5)) + scale_x_discrete(position = "top") +
  theme(axis.text.y = element_text(hjust = 1, size = 8)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 40, hjust = 0, size = 8)) +
  # theme(axis.text.x = element_blank()) +
  # theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  labs(alpha = "Adjusted P value", size = "NES abs. value", color = "Enriched samples") +
  theme(plot.margin = margin(0.2, 2, 0.2, 0.2, "cm"))
#theme(text = element_text(face = "bold"))
pdf(file = "plots/timecourse/timecourse-d7911-dotplots-hallmark-EuploidGATA1s-vs-EuploidwtGATA1-d7.pdf")  
p <- p + guides(fill = FALSE) + theme(legend.position = "right", legend.box = "vertical")
p <- p + guides(color = guide_legend(order = 1), size  = guide_legend(order = 2), alpha = guide_legend(order = 2))
print(p)
dev.off()


### EuploidGATAsD11 - EuploidwtGATA1D11 as baseline reference ####

# list dge sample cell type files
dge_files <- list.files(path = "results/timecourse/EuploidGATA1s_vs_EuploidwtGATA1/D11", 
                        pattern = "timecourse_d7911_d11_")

# run gsea on each dge file
for (dge_file in dge_files) {
  # load dge results
  dge_results <- readr::read_tsv(file.path("results/timecourse/EuploidGATA1s_vs_EuploidwtGATA1/D11", 
                                           dge_file))
  
  # pre-rank lists of fold change by sorting
  lfc_rank <- dge_results %>% 
    dplyr::pull(avg_log2FC, name = gene_symbol)
  lfc_rank <- sort(lfc_rank, decreasing = TRUE)
  
  # run gsea
  # gsea_results <- fgsea::fgsea(pathways = gene_set, stats = lfc_rank) %>%
  gsea_results <- fgsea::fgsea(pathways = gene_set, stats = lfc_rank, 
                               nPermSimple = 100000) %>% 
    tibble::as_tibble() %>% 
    dplyr::filter(padj <= 0.05) %>% 
    dplyr::arrange(dplyr::desc(NES))
  
  # write gsea results to file
  output_file <- paste0(strsplit(dge_file, 
                                 split = "_dge")[[1]][1], "_gsea.tsv.gz")
  gsea_results %>%  readr::write_tsv(file.path("results/timecourse/EuploidGATA1s_vs_EuploidwtGATA1/D11", 
                                               "GSEA", output_file))
}

# list gsea sample cell type files
gsea_dir <- file.path("results/timecourse/EuploidGATA1s_vs_EuploidwtGATA1/D11", "GSEA")
gsea_files <- list.files(path = gsea_dir, pattern = "_gsea.tsv.gz")

# initialize list for gsea results
gsea_results_list <-list()

for (gsea_file in gsea_files) {
  # load dge results
  contrast <- strsplit(gsea_file, split = "_gsea")[[1]][1] 
  contrast <- gsub("timecourse_d7911_d11_", "", contrast)
  gsea_results <- readr::read_tsv(file.path(gsea_dir, gsea_file)) %>% 
    dplyr::select(pathway, padj, NES) %>% 
    dplyr::rename(GeneSet = pathway) %>% 
    dplyr::mutate(GeneSet = gsub("HALLMARK_", "", GeneSet)) %>% 
    dplyr::mutate(GeneSet = gsub("_", " ", GeneSet)) %>%
    dplyr::mutate(Group = contrast) %>% 
    dplyr::mutate(Group = case_when(Group == "HPCs_EuploidGATA1sD11" ~ "HPCs",
                                    Group == "HPC_-_MK_bias_2_EuploidGATA1sD11" ~ "HPC - MK bias 2",
                                    Group == "HPC_-_MK_bias_1_EuploidGATA1sD11" ~ "HPC - MK bias 1",
                                    Group == "MK_EuploidGATA1sD11" ~ "MK",
                                    Group == "HPC_-_Ery_bias_EuploidGATA1sD11" ~ "HPC - Ery bias",
                                    Group == "Ery_EuploidGATA1sD11" ~ "Ery",
                                    Group == "Myeloid_EuploidGATA1sD11" ~ "Myeloid")) %>%
    dplyr::mutate(Enriched.group = case_when(NES > 0 ~ "1", NES < 0 ~ "2"))
  gsea_results_list <- append(gsea_results_list, list(gsea_results))
}

# merge gsea results
gsea_results_merged <- bind_rows(gsea_results_list)
gsea_results_merged$Group <- factor(gsea_results_merged$Group, 
                                    levels = c("HPCs", "HPC - MK bias 2", "HPC - MK bias 1", 
                                               "MK", "HPC - Ery bias", "Ery", "Myeloid"))

# dotplot
p <- ggplot(gsea_results_merged, aes(Group, GeneSet)) + 
  geom_point(aes(fill = Enriched.group, alpha = padj, shape =21, size = abs(NES))) + 
  geom_point(aes(color = Enriched.group, shape =21, size = abs(NES)))  +
  scale_shape_identity() 

p <- p + theme_bw() + scale_color_manual(values = c( "1" = "forestgreen", "2" = "royalblue"), 
                                         labels = c("1" = "Euploid/GATA1s", "2" = "Euploid/wtGATA1")) +
  scale_fill_manual(values = c("forestgreen", "royalblue")) + scale_alpha_continuous(range = c(1,0.1)) +
  scale_size(range = c(1,5)) + scale_x_discrete(position = "top") +
  theme(axis.text.y = element_text(hjust = 1, size = 8)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 40, hjust = 0, size = 8)) +
  # theme(axis.text.x = element_blank()) +
  # theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  labs(alpha = "Adjusted P value", size = "NES abs. value", color = "Enriched samples") +
  theme(plot.margin = margin(0.2, 2, 0.2, 0.2, "cm"))
#theme(text = element_text(face = "bold"))
pdf(file = "plots/timecourse/timecourse-d7911-dotplots-hallmark-EuploidGATA1s-vs-EuploidwtGATA1-d11.pdf")  
p <- p + guides(fill = FALSE) + theme(legend.position = "right", legend.box = "vertical")
p <- p + guides(color = guide_legend(order = 1), size  = guide_legend(order = 2), alpha = guide_legend(order = 2))
print(p)
dev.off()


### T21wtGATAD7 - EuploidwtGATA1D7 as baseline reference ####

# list dge sample cell type files
dge_files <- list.files(path = "results/timecourse/T21wtGATA1_vs_EuploidwtGATA1/D7", 
                        pattern = "timecourse_d7911_d7_")

# run gsea on each dge file
for (dge_file in dge_files) {
  # load dge results
  dge_results <- readr::read_tsv(file.path("results/timecourse/T21wtGATA1_vs_EuploidwtGATA1/D7", 
                                           dge_file))
  
  # pre-rank lists of fold change by sorting
  lfc_rank <- dge_results %>% 
    dplyr::pull(avg_log2FC, name = gene_symbol)
  lfc_rank <- sort(lfc_rank, decreasing = TRUE)
  
  # run gsea
  # gsea_results <- fgsea::fgsea(pathways = gene_set, stats = lfc_rank) %>%
  gsea_results <- fgsea::fgsea(pathways = gene_set, stats = lfc_rank, 
                               nPermSimple = 100000) %>% 
    tibble::as_tibble() %>% 
    dplyr::filter(padj <= 0.05) %>% 
    dplyr::arrange(dplyr::desc(NES))
  
  # write gsea results to file
  output_file <- paste0(strsplit(dge_file, 
                                 split = "_dge")[[1]][1], "_gsea.tsv.gz")
  gsea_results %>%  readr::write_tsv(file.path("results/timecourse/T21wtGATA1_vs_EuploidwtGATA1/D7", 
                                               "GSEA", output_file))
}

# list gsea sample cell type files
gsea_dir <- file.path("results/timecourse/T21wtGATA1_vs_EuploidwtGATA1/D7", "GSEA")
gsea_files <- list.files(path = gsea_dir, pattern = "_gsea.tsv.gz")

# initialize list for gsea results
gsea_results_list <-list()

for (gsea_file in gsea_files) {
  # load dge results
  contrast <- strsplit(gsea_file, split = "_gsea")[[1]][1] 
  contrast <- gsub("timecourse_d7911_d7_", "", contrast)
  gsea_results <- readr::read_tsv(file.path(gsea_dir, gsea_file)) %>% 
    dplyr::select(pathway, padj, NES) %>% 
    dplyr::rename(GeneSet = pathway) %>% 
    dplyr::mutate(GeneSet = gsub("HALLMARK_", "", GeneSet)) %>% 
    dplyr::mutate(GeneSet = gsub("_", " ", GeneSet)) %>%
    dplyr::mutate(Group = contrast) %>% 
    dplyr::mutate(Group = case_when(Group == "HPCs_T21wtGATA1D7" ~ "HPCs",
                                    Group == "HPC_-_MK_bias_1_T21wtGATA1D7" ~ "HPC - MK bias 1",
                                    Group == "MK_T21wtGATA1D7" ~ "MK")) %>%
    dplyr::mutate(Enriched.group = case_when(NES > 0 ~ "1", NES < 0 ~ "2"))
  gsea_results_list <- append(gsea_results_list, list(gsea_results))
}

# merge gsea results
gsea_results_merged <- bind_rows(gsea_results_list)
gsea_results_merged$Group <- factor(gsea_results_merged$Group, 
                                    levels = c("HPCs", "HPC - MK bias 1", "MK"))

# dotplot
p <- ggplot(gsea_results_merged, aes(Group, GeneSet)) + 
  geom_point(aes(fill = Enriched.group, alpha = padj, shape =21, size = abs(NES))) + 
  geom_point(aes(color = Enriched.group, shape =21, size = abs(NES)))  +
  scale_shape_identity() 

p <- p + theme_bw() + scale_color_manual(values = c( "1" = "forestgreen", "2" = "royalblue"), 
                                         labels = c("1" = "T21/wtGATA1", "2" = "Euploid/wtGATA1")) +
  scale_fill_manual(values = c("forestgreen", "royalblue")) + scale_alpha_continuous(range = c(1,0.1)) +
  scale_size(range = c(1,5)) + scale_x_discrete(position = "top") +
  theme(axis.text.y = element_text(hjust = 1, size = 8)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 40, hjust = 0, size = 8)) +
  # theme(axis.text.x = element_blank()) +
  # theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  labs(alpha = "Adjusted P value", size = "NES abs. value", color = "Enriched samples") +
  theme(plot.margin = margin(0.2, 2, 0.2, 0.2, "cm"))
#theme(text = element_text(face = "bold"))
pdf(file = "plots/timecourse/timecourse-d7911-dotplots-hallmark-T21wtGATA1_vs_EuploidwtGATA1-d7.pdf")  
p <- p + guides(fill = FALSE) + theme(legend.position = "right", legend.box = "vertical")
p <- p + guides(color = guide_legend(order = 1), size  = guide_legend(order = 2), alpha = guide_legend(order = 2))
print(p)
dev.off()


### T21wtGATAD11 - EuploidwtGATA1D11 as baseline reference ####

# list dge sample cell type files
dge_files <- list.files(path = "results/timecourse/T21wtGATA1_vs_EuploidwtGATA1/D11", 
                        pattern = "timecourse_d7911_d11_")

# run gsea on each dge file
for (dge_file in dge_files) {
  # load dge results
  dge_results <- readr::read_tsv(file.path("results/timecourse/T21wtGATA1_vs_EuploidwtGATA1/D11", 
                                           dge_file))
  
  # pre-rank lists of fold change by sorting
  lfc_rank <- dge_results %>% 
    dplyr::pull(avg_log2FC, name = gene_symbol)
  lfc_rank <- sort(lfc_rank, decreasing = TRUE)
  
  # run gsea
  # gsea_results <- fgsea::fgsea(pathways = gene_set, stats = lfc_rank) %>%
  gsea_results <- fgsea::fgsea(pathways = gene_set, stats = lfc_rank, 
                               nPermSimple = 100000) %>% 
    tibble::as_tibble() %>% 
    dplyr::filter(padj <= 0.05) %>% 
    dplyr::arrange(dplyr::desc(NES))
  
  # write gsea results to file
  output_file <- paste0(strsplit(dge_file, 
                                 split = "_dge")[[1]][1], "_gsea.tsv.gz")
  gsea_results %>%  readr::write_tsv(file.path("results/timecourse/T21wtGATA1_vs_EuploidwtGATA1/D11", 
                                               "GSEA", output_file))
}

# list gsea sample cell type files
gsea_dir <- file.path("results/timecourse/T21wtGATA1_vs_EuploidwtGATA1/D11", "GSEA")
gsea_files <- list.files(path = gsea_dir, pattern = "_gsea.tsv.gz")

# initialize list for gsea results
gsea_results_list <-list()

for (gsea_file in gsea_files) {
  # load dge results
  contrast <- strsplit(gsea_file, split = "_gsea")[[1]][1] 
  contrast <- gsub("timecourse_d7911_d11_", "", contrast)
  gsea_results <- readr::read_tsv(file.path(gsea_dir, gsea_file)) %>% 
    dplyr::select(pathway, padj, NES) %>% 
    dplyr::rename(GeneSet = pathway) %>% 
    dplyr::mutate(GeneSet = gsub("HALLMARK_", "", GeneSet)) %>% 
    dplyr::mutate(GeneSet = gsub("_", " ", GeneSet)) %>%
    dplyr::mutate(Group = contrast) %>% 
    dplyr::mutate(Group = case_when(Group == "HPCs_T21wtGATA1D11" ~ "HPCs",
                                    Group == "HPC_-_MK_bias_2_T21wtGATA1D11" ~ "HPC - MK bias 2",
                                    Group == "HPC_-_MK_bias_1_T21wtGATA1D11" ~ "HPC - MK bias 1",
                                    Group == "MK_T21wtGATA1D11" ~ "MK",
                                    Group == "Ery_T21wtGATA1D11" ~ "Ery",
                                    Group == "Myeloid_T21wtGATA1D11" ~ "Myeloid")) %>%
    dplyr::mutate(Enriched.group = case_when(NES > 0 ~ "1", NES < 0 ~ "2"))
  gsea_results_list <- append(gsea_results_list, list(gsea_results))
}

# merge gsea results
gsea_results_merged <- bind_rows(gsea_results_list)
gsea_results_merged$Group <- factor(gsea_results_merged$Group, 
                                    levels = c("HPCs", "HPC - MK bias 2", "HPC - MK bias 1", 
                                               "MK", "Ery", "Myeloid"))

# dotplot
p <- ggplot(gsea_results_merged, aes(Group, GeneSet)) + 
  geom_point(aes(fill = Enriched.group, alpha = padj, shape =21, size = abs(NES))) + 
  geom_point(aes(color = Enriched.group, shape =21, size = abs(NES)))  +
  scale_shape_identity() 

p <- p + theme_bw() + scale_color_manual(values = c( "1" = "forestgreen", "2" = "royalblue"), 
                                         labels = c("1" = "T21/wtGATA1", "2" = "Euploid/wtGATA1")) +
  scale_fill_manual(values = c("forestgreen", "royalblue")) + scale_alpha_continuous(range = c(1,0.1)) +
  scale_size(range = c(1,5)) + scale_x_discrete(position = "top") +
  theme(axis.text.y = element_text(hjust = 1, size = 8)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 40, hjust = 0, size = 8)) +
  # theme(axis.text.x = element_blank()) +
  # theme(axis.ticks.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  labs(alpha = "Adjusted P value", size = "NES abs. value", color = "Enriched samples") +
  theme(plot.margin = margin(0.2, 2, 0.2, 0.2, "cm"))
#theme(text = element_text(face = "bold"))
pdf(file = "plots/timecourse/timecourse-d7911-dotplots-hallmark-T21wtGATA1-vs-EuploidwtGATA1-d11.pdf")  
p <- p + guides(fill = FALSE) + theme(legend.position = "right", legend.box = "vertical")
p <- p + guides(color = guide_legend(order = 1), size  = guide_legend(order = 2), alpha = guide_legend(order = 2))
print(p)
dev.off()

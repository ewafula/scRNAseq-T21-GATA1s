# Prepare data from D7/D9/D11 timecourse Seurat object for MAST
# differential gene expression. Asking the question:
# What genes could promote a stalled phenotype unique to the T21/GATA1s 
# condition in early erythroid, megakaryocyte, and myeloid cell annotations.
# Reference - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4362246/pdf/JCI75714.pdf 


# Eric Wafula
# 2024

# Load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(MAST))

# set seed for reproducibility
set.seed(123)

# establish directories
RDS_dir <- file.path("..", "seurat_objects")
results_dir <- file.path("..", "results", "timecourse")

# load integrated and  annotated Seurat RDS Object
integrated <- readRDS(file.path(RDS_dir, "timecourse-d7911-renamed.RDS"))

# update column metadata
integrated@meta.data <- integrated@meta.data %>% 
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
                                 grepl("D11", orig.ident) ~ "D11"),
                lineage = case_when(cell_type == "HPCs"  ~ "HPCs",
                                    cell_type == "HPC - MK bias 2" ~ "MK",
                                    cell_type == "HPC - MK bias 1" ~ "MK",
                                    cell_type == "MK" ~ "MK",
                                    cell_type == "HPC - Ery bias" ~ "Ery",
                                    cell_type == "Ery" ~ "Ery",
                                    cell_type == "Myeloid" ~ "Myeloid"))
integrated@meta.data$sample <- integrated@meta.data$orig.ident
integrated@meta.data$barcode <- rownames(integrated@meta.data)

# add scVelo velocity latent time 
latent_time <- readr::read_csv(file.path(results_dir, 
                               "timecourse-d7911-renamed-latent-time.csv.gz"))
barcodes = rownames(integrated@meta.data)
latent_time = latent_time[match(barcodes, latent_time$barcode), ] #match the order of seurat barcode
integrated@meta.data$latent_time = latent_time$latent_time #put laten time into metadata

# back transform SCT assay
integrated <- Seurat::PrepSCTFindMarkers(integrated, assay = "SCT")

# convert Seurat object to SingleCellExperiment object keeping the SCT assay
Seurat::DefaultAssay(integrated) <- "SCT" 
sce <- as.SingleCellExperiment(
  Seurat::DietSeurat(integrated, assays = "SCT", graphs = c("pca", "umap")))

# convert  SingleCellExperiment object to MAST SingleCellAssay object
sca <- MAST::SceToSingleCellAssay(sce, class = "SingleCellAssay")

# subset for HVGs that were used in integration
hvg <-
  readr::read_lines(file.path(results_dir,
                              "timecourse-d7911-umap-renamed-variable-genes.txt")) 
sca <- sca[hvg, ]
saveRDS(sca, 
        file.path(RDS_dir, "timecourse-d7911-renamed-sca.RDS"))

###### DE for the "HPC - Ery bias"
# subset sca object
hpc_ery_bias_sca <- sca[, sca@colData$cell_type == "HPC - Ery bias"]
saveRDS(hpc_ery_bias_sca,
        file.path(RDS_dir, "timecourse-d7911-renamed-HPC-Ery-bias-sca.RDS"))

# fit a hurdle model, modeling the condition, genotype, and their interaction and
# latent_time fator, thus adjusting for pseudotemporal cell ordering along the lineage
cond <- factor(x = colData(hpc_ery_bias_sca)$condition)
cond <- relevel(x = cond, ref = "Euploid")
colData(hpc_ery_bias_sca)$condition <- cond
cond <- factor(x = colData(hpc_ery_bias_sca)$genotype)
cond <- relevel(x = cond, ref = "wtGATA1")
colData(hpc_ery_bias_sca)$genotype <- cond
fmla <- as.formula("~ condition*genotype + latent_time")
# options(mc.cores = 10)
zlmCond <- MAST::zlm(formula = fmla, sca = hpc_ery_bias_sca, parallel=TRUE)
saveRDS(zlmCond,
        file.path(RDS_dir, "timecourse-d7911-renamed-mast-zlmCond-HPC-Ery-bias.RDS"))

# summarize
summaryCond <- MAST::summary(object = zlmCond,
                             doLRT = 'conditionT21:genotypeGATA1s',
                             parallel=TRUE)
summaryH <- summaryCond$datatable %>%
  dplyr::filter(component == "H") %>%
  dplyr::select(primerid, `Pr(>Chisq)`)
summaryFC <- summaryCond$datatable %>%
  dplyr::filter(component == "logFC") %>%
  dplyr::select(primerid, contrast, coef) %>%
  tidyr::pivot_wider(names_from = contrast, values_from = coef)
summaryDt <- dplyr::left_join(summaryH, summaryFC, by=join_by(primerid)) %>%
  dplyr::arrange(`Pr(>Chisq)`)
summaryDt[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
summaryDt %>%
  readr::write_tsv(file.path(results_dir,
                             "timecourse-d7911-renamed-mast-zlmCond-HPC-Ery-bias-DEGs.tsv.gz"))

###### DE for the "HPC - MK bias 1/2"
# subset sca object
hpc_mk_bias_sca <- sca[, sca@colData$cell_type %in% c("HPC - MK bias 2", "HPC - MK bias 1")]
saveRDS(hpc_mk_bias_sca,
        file.path(RDS_dir, "timecourse-d7911-renamed-HPC-MK-bias-sca.RDS"))

# fit a hurdle model, modeling the condition, genotype, and their interaction and
# latent_time fator, thus adjusting for pseudotemporal cell ordering along the lineage
cond <- factor(x = colData(hpc_mk_bias_sca)$condition)
cond <- relevel(x = cond, ref = "Euploid")
colData(hpc_mk_bias_sca)$condition <- cond
cond <- factor(x = colData(hpc_mk_bias_sca)$genotype)
cond <- relevel(x = cond, ref = "wtGATA1")
colData(hpc_mk_bias_sca)$genotype <- cond
fmla <- as.formula("~ condition*genotype + latent_time")
# options(mc.cores = 10)
zlmCond <- MAST::zlm(formula = fmla, sca = hpc_mk_bias_sca, parallel=TRUE)
saveRDS(zlmCond,
        file.path(RDS_dir, "timecourse-d7911-renamed-mast-zlmCond-HPC-MK-bias.RDS"))

# summarize
summaryCond <- MAST::summary(object = zlmCond,
                             doLRT = 'conditionT21:genotypeGATA1s',
                             parallel=TRUE)
summaryH <- summaryCond$datatable %>%
  dplyr::filter(component == "H") %>%
  dplyr::select(primerid, `Pr(>Chisq)`)
summaryFC <- summaryCond$datatable %>%
  dplyr::filter(component == "logFC") %>%
  dplyr::select(primerid, contrast, coef) %>%
  tidyr::pivot_wider(names_from = contrast, values_from = coef)
summaryDt <- dplyr::left_join(summaryH, summaryFC, by=join_by(primerid)) %>%
  dplyr::arrange(`Pr(>Chisq)`)
summaryDt[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
summaryDt %>%
  readr::write_tsv(file.path(results_dir,
                             "timecourse-d7911-renamed-mast-zlmCond-HPC-MK-bias-DEGs.tsv.gz"))

sca <- readRDS(file.path(RDS_dir, "timecourse-d7911-renamed-sca.RDS"))

###### DE for the "HPCs"
# subset sca object
hpcs_sca <- sca[, sca@colData$cell_type %in% c("HPCs")]
saveRDS(hpcs_sca, 
        file.path(RDS_dir, "timecourse-d7911-renamed-HPCs-sca.RDS"))

# fit a hurdle model, modeling the condition, genotype, and their interaction and
# latent_time fator, thus adjusting for pseudotemporal cell ordering along the lineage
cond <- factor(x = colData(hpcs_sca)$condition)
cond <- relevel(x = cond, ref = "Euploid")
colData(hpcs_sca)$condition <- cond
cond <- factor(x = colData(hpcs_sca)$genotype)
cond <- relevel(x = cond, ref = "wtGATA1")
colData(hpcs_sca)$genotype <- cond
fmla <- as.formula("~ condition*genotype + latent_time")
# options(mc.cores = 10)
zlmCond <- MAST::zlm(formula = fmla, sca = hpcs_sca, parallel=TRUE)
saveRDS(zlmCond, 
        file.path(RDS_dir, "timecourse-d7911-renamed-mast-zlmCond-HPCs.RDS"))

# summarize
summaryCond <- MAST::summary(object = zlmCond, 
                             doLRT = 'conditionT21:genotypeGATA1s', 
                             parallel=TRUE)
summaryH <- summaryCond$datatable %>% 
  dplyr::filter(component == "H") %>% 
  dplyr::select(primerid, `Pr(>Chisq)`)
summaryFC <- summaryCond$datatable %>% 
  dplyr::filter(component == "logFC") %>% 
  dplyr::select(primerid, contrast, coef) %>% 
  tidyr::pivot_wider(names_from = contrast, values_from = coef)
summaryDt <- dplyr::left_join(summaryH, summaryFC, by=join_by(primerid)) %>% 
  dplyr::arrange(`Pr(>Chisq)`)
summaryDt[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
summaryDt %>%
  readr::write_tsv(file.path(results_dir,
                             "timecourse-d7911-renamed-mast-zlmCond-HPCs-DEGs.tsv.gz"))


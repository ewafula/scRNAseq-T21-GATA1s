# Single-cell transcriptomics reveal individual and cooperative effects of Trisomy 21 and GATA1s on hematopoiesis


This repository contains the scripts used for the analysis shown in [Single-cell transcriptomics reveal individual and cooperative effects of Trisomy 21 and GATA1s on hematopoiesis]().

## Summary
Trisomy 21 (T21) is associated with baseline erythrocytosis, thrombocytopenia, neutrophilia, transient abnormal myelopoiesis (TAM), and myeloid leukemia of Down syndrome (ML-DS). TAM and ML-DS blasts harbor mutations in GATA1, resulting in the exclusive expression of the truncated isoform GATA1s. Germline GATA1s mutations in individuals without T21 cause congenital cytopenias, typically without a leukemic predisposition. To dissect the developmental effects of T21 and GATA1s, we used a combination of isogenic human induced pluripotent stem cells (iPSCs), primary human fetal and neonatal cells, and single-cell transcriptomics to interrogate hematopoietic progenitors differing only by chromosome 21 and/or GATA1 status. Both T21 and GATA1s induced early lineage skewing, and trajectory analysis revealed that GATA1s altered the temporal regulation of lineage-specific transcriptional programs, disrupting cell proliferation and maturation irrespective of chromosomal context. These studies uncovered unexpected heterogeneity and lineage priming in early, multipotent hematopoietic progenitors, and identified transcriptional and functional maturation blocks linked to GATA1s.



## Analysis workflow directory structure
This illustration shows the directory structure of the analysis workflow, which can be used to reproduce the study's analysis results in the publication.

Please note that the software packages and code implementations may change when using this archived repository. You may need to update the code to be compatible with the latest versions of the software packages.

The 10x scRNA-Seq sequencing data, external reference Rdata, and the analysis Seurat objects have not been uploaded to the repository due to file size limitations. 
```
.
|-- 01-run-soupx.R
|-- 02-run-doubletfinder.R
|-- 03-d7-and-primary-integrate-samples.R
|-- 03-d7-and-primary-reference-mapping.R
|-- 03-d7911-integrate-samples.R
|-- 04-d7-and-primary-annotate-clusters.R
|-- 04-d7911-annotate-clusters.R
|-- 05-d7911-gene-set-enrichment.R
|-- README.md
|-- data
|   |-- primary
|   |   |-- HTS_SO48_01_FL21
|   |   |-- HTS_SO48_02_FL33
|   |   |-- HTS_SO48_03_FL38
|   |   |-- HTS_SO48_04_FL50
|   |   |-- HTS_SO48_05_FL54
|   |   |-- HTS_SO48_06_TMD145
|   |   |-- HTS_SO48_07_TMD160
|   |   |-- HTS_SO92_13_FL32_T21_CD34plus
|   |   `-- data.info
|   `-- timecourse
|       |-- DS145_EM9_D11
|       |-- DS145_EM9_D7
|       |-- DS145_EM9_D9
|       |-- TMD145_EM331434_D11
|       |-- TMD145_EM331434_D7
|       |-- TMD145_EM331434_D9
|       |-- TMD145_EM3319_D11
|       |-- TMD145_EM3319_D7
|       |-- TMD145_EM3319_D9
|       |-- TMD145_SCM71_D11
|       |-- TMD145_SCM71_D7
|       |-- TMD145_SCM71_D9
|       `-- data.info
|-- metadata
|   |-- query_primary_plus_d7.tsv.gz
|   |-- reference.RDS.tsv.gz
|   |-- timecourse-d7-and-primary-renamed.tsv.gz
|   `-- ttimecourse-d7911-renamed.tsv.gz
|-- plots
|   |-- doubletfinder
|   |   |-- DS145_EM9_D11_pk_sweep_plot.pdf
|   |   |-- DS145_EM9_D11_umap_doublets.pdf
|   |   |-- DS145_EM9_D7_pk_sweep_plot.pdf
|   |   |-- DS145_EM9_D7_umap_doublets.pdf
|   |   |-- DS145_EM9_D9_pk_sweep_plot.pdf
|   |   |-- DS145_EM9_D9_umap_doublets.pdf
|   |   |-- HTS_SO48_01_FL21_pk_sweep_plot.pdf
|   |   |-- HTS_SO48_01_FL21_umap_doublets.pdf
|   |   |-- HTS_SO48_02_FL33_pk_sweep_plot.pdf
|   |   |-- HTS_SO48_02_FL33_umap_doublets.pdf
|   |   |-- HTS_SO48_03_FL38_pk_sweep_plot.pdf
|   |   |-- HTS_SO48_03_FL38_umap_doublets.pdf
|   |   |-- HTS_SO48_04_FL50_pk_sweep_plot.pdf
|   |   |-- HTS_SO48_04_FL50_umap_doublets.pdf
|   |   |-- HTS_SO48_05_FL54_pk_sweep_plot.pdf
|   |   |-- HTS_SO48_05_FL54_umap_doublets.pdf
|   |   |-- HTS_SO48_06_TMD145_pk_sweep_plot.pdf
|   |   |-- HTS_SO48_06_TMD145_umap_doublets.pdf
|   |   |-- HTS_SO48_07_TMD160_pk_sweep_plot.pdf
|   |   |-- HTS_SO48_07_TMD160_umap_doublets.pdf
|   |   |-- TMD145_EM331434_D11_pk_sweep_plot.pdf
|   |   |-- TMD145_EM331434_D11_umap_doublets.pdf
|   |   |-- TMD145_EM331434_D7_pk_sweep_plot.pdf
|   |   |-- TMD145_EM331434_D7_umap_doublets.pdf
|   |   |-- TMD145_EM331434_D9_pk_sweep_plot.pdf
|   |   |-- TMD145_EM331434_D9_umap_doublets.pdf
|   |   |-- TMD145_EM3319_D11_pk_sweep_plot.pdf
|   |   |-- TMD145_EM3319_D11_umap_doublets.pdf
|   |   |-- TMD145_EM3319_D7_pk_sweep_plot.pdf
|   |   |-- TMD145_EM3319_D7_umap_doublets.pdf
|   |   |-- TMD145_EM3319_D9_pk_sweep_plot.pdf
|   |   |-- TMD145_EM3319_D9_umap_doublets.pdf
|   |   |-- TMD145_SCM71_D11_pk_sweep_plot.pdf
|   |   |-- TMD145_SCM71_D11_umap_doublets.pdf
|   |   |-- TMD145_SCM71_D7_pk_sweep_plot.pdf
|   |   |-- TMD145_SCM71_D7_umap_doublets.pdf
|   |   |-- TMD145_SCM71_D9_pk_sweep_plot.pdf
|   |   `-- TMD145_SCM71_D9_umap_doublets.pdf
|   |-- primary
|   |   |-- timecourse-d7-and-primary-cluster-annotation.pdf
|   |   |-- timecourse-d7-and-primary-elbowplot.pdf
|   |   |-- timecourse-d7-and-primary-ordering1-ery.pdf
|   |   |-- timecourse-d7-and-primary-ordering1-mk.pdf
|   |   |-- timecourse-d7-and-primary-ordering1-myeloid.pdf
|   |   |-- timecourse-d7-and-primary-ordering2-ery-violin.pdf
|   |   |-- timecourse-d7-and-primary-ordering2-ery.pdf
|   |   |-- timecourse-d7-and-primary-ordering2-mk-violin.pdf
|   |   |-- timecourse-d7-and-primary-ordering2-mk.pdf
|   |   |-- timecourse-d7-and-primary-ordering2-myeloid-violin.pdf
|   |   |-- timecourse-d7-and-primary-ordering2-myeloid.pdf
|   |   |-- timecourse-d7-and-primary-scorecard1-ery.pdf
|   |   |-- timecourse-d7-and-primary-scorecard1-mk.pdf
|   |   |-- timecourse-d7-and-primary-scorecard2-ery.pdf
|   |   |-- timecourse-d7-and-primary-scorecard2-mk.pdf
|   |   |-- timecourse-d7-and-primary-umap-marker-dotplots.pdf
|   |   |-- timecourse-d7-and-primary-umap-renamed-proportions.pdf
|   |   |-- timecourse-d7-and-primary-umap-renamed-split.pdf
|   |   |-- timecourse-d7-and-primary-umap-renamed-v2.pdf
|   |   |-- timecourse-d7-and-primary-umap-renamed.pdf
|   |   |-- timecourse-d7-and-primary-umap-sample_type.pdf
|   |   |-- timecourse-d7-and-primary-umap.pdf
|   |   `-- timecourse-d7-and-primary-umap_bigfont.pdf
|   |-- reference
|   |   |-- query_umap_predicted_cell_types_primary_plus_d7-proportions.pdf
|   |   |-- query_umap_predicted_cell_types_primary_plus_d7.pdf
|   |   |-- query_umap_sample_types_primary_plus_d7.pdf
|   |   |-- query_umap_samples_primary_plus_d7-ordering1-ery.pdf
|   |   |-- query_umap_samples_primary_plus_d7-ordering1-mk.pdf
|   |   |-- query_umap_samples_primary_plus_d7-ordering2-ery.pdf
|   |   |-- query_umap_samples_primary_plus_d7-ordering2-mk.pdf
|   |   |-- query_umap_samples_primary_plus_d7_scorecard1_ery.pdf
|   |   |-- query_umap_samples_primary_plus_d7_scorecard1_mk.pdf
|   |   |-- query_umap_samples_primary_plus_d7_scorecard2_ery.pdf
|   |   |-- query_umap_samples_primary_plus_d7_scorecard2_mk.pdf
|   |   |-- query_umap_samples_split_primary_plus_d7.pdf
|   |   |-- reference_query_umap_sample_types_primary_plus_d7.pdf
|   |   |-- reference_query_umap_samples_primary_plus_d7.pdf
|   |   |-- reference_umap_cell_types.pdf
|   |   |-- reference_umap_cell_types_v2.pdf
|   |   |-- reference_umap_samples.pdf
|   |   `-- reference_umap_samples_split.pdf
|   |-- soupx
|   |   |-- DS145_EM9_D11_soupx_contaminant.pdf
|   |   |-- DS145_EM9_D7_soupx_contaminant.pdf
|   |   |-- DS145_EM9_D9_soupx_contaminant.pdf
|   |   |-- HTS_SO48_01_FL21_soupx_contaminant.pdf
|   |   |-- HTS_SO48_02_FL33_soupx_contaminant.pdf
|   |   |-- HTS_SO48_03_FL38_soupx_contaminant.pdf
|   |   |-- HTS_SO48_04_FL50_soupx_contaminant.pdf
|   |   |-- HTS_SO48_05_FL54_soupx_contaminant.pdf
|   |   |-- HTS_SO48_06_TMD145_soupx_contaminant.pdf
|   |   |-- HTS_SO48_07_TMD160_soupx_contaminant.pdf
|   |   |-- TMD145_EM331434_D11_soupx_contaminant.pdf
|   |   |-- TMD145_EM331434_D7_soupx_contaminant.pdf
|   |   |-- TMD145_EM331434_D9_soupx_contaminant.pdf
|   |   |-- TMD145_EM3319_D11_soupx_contaminant.pdf
|   |   |-- TMD145_EM3319_D7_soupx_contaminant.pdf
|   |   |-- TMD145_EM3319_D9_soupx_contaminant.pdf
|   |   |-- TMD145_SCM71_D11_soupx_contaminant.pdf
|   |   |-- TMD145_SCM71_D7_soupx_contaminant.pdf
|   |   `-- TMD145_SCM71_D9_soupx_contaminant.pdf
|   `-- timecourse
|       |-- timecourse-d7911-dotplots-hallmark-EuploidGATA1s-vs-EuploidwtGATA1-d11.pdf
|       |-- timecourse-d7911-dotplots-hallmark-EuploidGATA1s-vs-EuploidwtGATA1-d7.pdf
|       |-- timecourse-d7911-dotplots-hallmark-T21GATA1s-vs-T21wtGATA1-d11.pdf
|       |-- timecourse-d7911-dotplots-hallmark-T21GATA1s-vs-T21wtGATA1-d7.pdf
|       |-- timecourse-d7911-dotplots-hallmark-T21wtGATA1-vs-EuploidwtGATA1-d11.pdf
|       |-- timecourse-d7911-dotplots-hallmark-T21wtGATA1_vs_EuploidwtGATA1-d7.pdf
|       |-- timecourse-d7911-elbowplot.pdf
|       |-- timecourse-d7911-erythroids-key-genes-umap.pdf
|       |-- timecourse-d7911-erythroids-umap.pdf
|       |-- timecourse-d7911-key-genes-doplot.pdf
|       |-- timecourse-d7911-megakaryocytes-key-genes-umap.pdf
|       |-- timecourse-d7911-megakaryocytes-umap.pdf
|       |-- timecourse-d7911-myeloids-key-genes-umap.pdf
|       |-- timecourse-d7911-myeloids-umap.pdf
|       |-- timecourse-d7911-ordering1-ery.pdf
|       |-- timecourse-d7911-ordering1-mk.pdf
|       |-- timecourse-d7911-ordering1-myeloid.pdf
|       |-- timecourse-d7911-ordering2-ery.pdf
|       |-- timecourse-d7911-ordering2-mk.pdf
|       |-- timecourse-d7911-ordering2-myeloid.pdf
|       |-- timecourse-d7911-ordering3-ery-violin.pdf
|       |-- timecourse-d7911-ordering3-ery.pdf
|       |-- timecourse-d7911-ordering3-hpc-violin.pdf
|       |-- timecourse-d7911-ordering3-hpcs.pdf
|       |-- timecourse-d7911-ordering3-mk-violin.pdf
|       |-- timecourse-d7911-ordering3-mk.pdf
|       |-- timecourse-d7911-ordering3-myeloid-violin.pdf
|       |-- timecourse-d7911-ordering3-myeloid.pdf
|       |-- timecourse-d7911-umap-renamed-marker-dotplots.pdf
|       |-- timecourse-d7911-umap-renamed-proportions-cell-cycle-phases.pdf
|       |-- timecourse-d7911-umap-renamed-proportions-d11.pdf
|       |-- timecourse-d7911-umap-renamed-proportions-d7.pdf
|       |-- timecourse-d7911-umap-renamed-proportions-d9.pdf
|       |-- timecourse-d7911-umap-renamed-split-v2.pdf
|       |-- timecourse-d7911-umap-renamed-split.pdf
|       |-- timecourse-d7911-umap-renamed.pdf
|       |-- timecourse-d7911-umap-renamed.v2.pdf
|       `-- timecourse-d7911-umap.pdf
|-- reference_data
|   |-- Rdata
|   |   |-- data.info
|   |   |-- seurat_object.Rdata
|   |   `-- seurat_object_rerun.Rdata
|   |-- data.info
|   |-- plots
|   |   |-- quality_by_origin.pdf
|   |   |-- tsne_clusterRX105.pdf
|   |   |-- tsne_originRX105.pdf
|   |   |-- umap_clusterRX105.pdf
|   |   `-- umap_originRX105.pdf
|   `-- seurat_analysis.R
|-- results
|   |-- doubletfinder
|   |   |-- DS145_EM9_D11_pk_values.txt
|   |   |-- DS145_EM9_D7_pk_values.txt
|   |   |-- DS145_EM9_D9_pk_values.txt
|   |   |-- HTS_SO48_01_FL21_pk_values.txt
|   |   |-- HTS_SO48_02_FL33_pk_values.txt
|   |   |-- HTS_SO48_03_FL38_pk_values.txt
|   |   |-- HTS_SO48_04_FL50_pk_values.txt
|   |   |-- HTS_SO48_05_FL54_pk_values.txt
|   |   |-- HTS_SO48_06_TMD145_pk_values.txt
|   |   |-- HTS_SO48_07_TMD160_pk_values.txt
|   |   |-- TMD145_EM331434_D11_pk_values.txt
|   |   |-- TMD145_EM331434_D7_pk_values.txt
|   |   |-- TMD145_EM331434_D9_pk_values.txt
|   |   |-- TMD145_EM3319_D11_pk_values.txt
|   |   |-- TMD145_EM3319_D7_pk_values.txt
|   |   |-- TMD145_EM3319_D9_pk_values.txt
|   |   |-- TMD145_SCM71_D11_pk_values.txt
|   |   |-- TMD145_SCM71_D7_pk_values.txt
|   |   `-- TMD145_SCM71_D9_pk_values.txt
|   |-- primary
|   |   `-- timecourse-d7-and-primary-umap-renamed-cell-counts.tsv.gz
|   |-- reference
|   |   `-- primary_plus_and_d7_query_sample_cell_type_counts.tsv
|   `-- timecourse
|       |-- EuploidGATA1s_vs_EuploidwtGATA1
|       |   |-- D11
|       |   |   |-- GSEA
|       |   |   |   |-- timecourse_d7911_d11_Ery_EuploidGATA1sD11_gsea.tsv.gz
|       |   |   |   |-- timecourse_d7911_d11_HPC_-_Ery_bias_EuploidGATA1sD11_gsea.tsv.gz
|       |   |   |   |-- timecourse_d7911_d11_HPC_-_MK_bias_1_EuploidGATA1sD11_gsea.tsv.gz
|       |   |   |   |-- timecourse_d7911_d11_HPC_-_MK_bias_2_EuploidGATA1sD11_gsea.tsv.gz
|       |   |   |   |-- timecourse_d7911_d11_HPCs_EuploidGATA1sD11_gsea.tsv.gz
|       |   |   |   |-- timecourse_d7911_d11_MK_EuploidGATA1sD11_gsea.tsv.gz
|       |   |   |   `-- timecourse_d7911_d11_Myeloid_EuploidGATA1sD11_gsea.tsv.gz
|       |   |   |-- timecourse_d7911_d11_Ery_EuploidGATA1sD11_dge.tsv.gz
|       |   |   |-- timecourse_d7911_d11_HPC_-_Ery_bias_EuploidGATA1sD11_dge.tsv.gz
|       |   |   |-- timecourse_d7911_d11_HPC_-_MK_bias_1_EuploidGATA1sD11_dge.tsv.gz
|       |   |   |-- timecourse_d7911_d11_HPC_-_MK_bias_2_EuploidGATA1sD11_dge.tsv.gz
|       |   |   |-- timecourse_d7911_d11_HPCs_EuploidGATA1sD11_dge.tsv.gz
|       |   |   |-- timecourse_d7911_d11_MK_EuploidGATA1sD11_dge.tsv.gz
|       |   |   `-- timecourse_d7911_d11_Myeloid_EuploidGATA1sD11_dge.tsv.gz
|       |   `-- D7
|       |       |-- GSEA
|       |       |   |-- timecourse_d7911_d7_Ery_EuploidGATA1sD7_gsea.tsv.gz
|       |       |   |-- timecourse_d7911_d7_HPC_-_Ery_bias_EuploidGATA1sD7_gsea.tsv.gz
|       |       |   |-- timecourse_d7911_d7_HPC_-_MK_bias_1_EuploidGATA1sD7_gsea.tsv.gz
|       |       |   |-- timecourse_d7911_d7_HPC_-_MK_bias_2_EuploidGATA1sD7_gsea.tsv.gz
|       |       |   |-- timecourse_d7911_d7_HPCs_EuploidGATA1sD7_gsea.tsv.gz
|       |       |   |-- timecourse_d7911_d7_MK_EuploidGATA1sD7_gsea.tsv.gz
|       |       |   `-- timecourse_d7911_d7_Myeloid_EuploidGATA1sD7_gsea.tsv.gz
|       |       |-- timecourse_d7911_d7_Ery_EuploidGATA1sD7_dge.tsv.gz
|       |       |-- timecourse_d7911_d7_HPC_-_Ery_bias_EuploidGATA1sD7_dge.tsv.gz
|       |       |-- timecourse_d7911_d7_HPC_-_MK_bias_1_EuploidGATA1sD7_dge.tsv.gz
|       |       |-- timecourse_d7911_d7_HPC_-_MK_bias_2_EuploidGATA1sD7_dge.tsv.gz
|       |       |-- timecourse_d7911_d7_HPCs_EuploidGATA1sD7_dge.tsv.gz
|       |       |-- timecourse_d7911_d7_MK_EuploidGATA1sD7_dge.tsv.gz
|       |       `-- timecourse_d7911_d7_Myeloid_EuploidGATA1sD7_dge.tsv.gz
|       |-- T21GATA1s_vs_T21wtGATA1
|       |   |-- D11
|       |   |   |-- GSEA
|       |   |   |   |-- timecourse_d7911_d11_Ery_T21GATA1sD11_gsea.tsv.gz
|       |   |   |   |-- timecourse_d7911_d11_HPC_-_Ery_bias_T21GATA1sD11_gsea.tsv.gz
|       |   |   |   |-- timecourse_d7911_d11_HPC_-_MK_bias_1_T21GATA1sD11_gsea.tsv.gz
|       |   |   |   |-- timecourse_d7911_d11_HPC_-_MK_bias_2_T21GATA1sD11_gsea.tsv.gz
|       |   |   |   `-- timecourse_d7911_d11_MK_T21GATA1sD11_gsea.tsv.gz
|       |   |   |-- timecourse_d7911_d11_Ery_T21GATA1sD11_dge.tsv.gz
|       |   |   |-- timecourse_d7911_d11_HPC_-_Ery_bias_T21GATA1sD11_dge.tsv.gz
|       |   |   |-- timecourse_d7911_d11_HPC_-_MK_bias_1_T21GATA1sD11_dge.tsv.gz
|       |   |   |-- timecourse_d7911_d11_HPC_-_MK_bias_2_T21GATA1sD11_dge.tsv.gz
|       |   |   |-- timecourse_d7911_d11_HPCs_T21GATA1sD11_dge.tsv.gz
|       |   |   |-- timecourse_d7911_d11_MK_T21GATA1sD11_dge.tsv.gz
|       |   |   `-- timecourse_d7911_d11_Myeloid_T21GATA1sD11_dge.tsv.gz
|       |   `-- D7
|       |       |-- GSEA
|       |       |   |-- timecourse_d7911_d7_HPC_-_MK_bias_1_T21GATA1sD7_gsea.tsv.gz
|       |       |   |-- timecourse_d7911_d7_HPC_-_MK_bias_2_T21GATA1sD7_gsea.tsv.gz
|       |       |   |-- timecourse_d7911_d7_HPCs_T21GATA1sD7_gsea.tsv.gz
|       |       |   |-- timecourse_d7911_d7_MK_T21GATA1sD7_gsea.tsv.gz
|       |       |   `-- timecourse_d7911_d7_Myeloid_T21GATA1sD7_gsea.tsv.gz
|       |       |-- timecourse_d7911_d7_Ery_T21GATA1sD7_dge.tsv.gz
|       |       |-- timecourse_d7911_d7_HPC_-_Ery_bias_T21GATA1sD7_dge.tsv.gz
|       |       |-- timecourse_d7911_d7_HPC_-_MK_bias_1_T21GATA1sD7_dge.tsv.gz
|       |       |-- timecourse_d7911_d7_HPC_-_MK_bias_2_T21GATA1sD7_dge.tsv.gz
|       |       |-- timecourse_d7911_d7_HPCs_T21GATA1sD7_dge.tsv.gz
|       |       |-- timecourse_d7911_d7_MK_T21GATA1sD7_dge.tsv.gz
|       |       `-- timecourse_d7911_d7_Myeloid_T21GATA1sD7_dge.tsv.gz
|       |-- T21wtGATA1_vs_EuploidwtGATA1
|       |   |-- D11
|       |   |   |-- GSEA
|       |   |   |-- timecourse_d7911_d11_Ery_T21wtGATA1D11_dge.tsv.gz
|       |   |   |-- timecourse_d7911_d11_HPC_-_Ery_bias_T21wtGATA1D11_dge.tsv.gz
|       |   |   |-- timecourse_d7911_d11_HPC_-_MK_bias_1_T21wtGATA1D11_dge.tsv.gz
|       |   |   |-- timecourse_d7911_d11_HPC_-_MK_bias_2_T21wtGATA1D11_dge.tsv.gz
|       |   |   |-- timecourse_d7911_d11_HPCs_T21wtGATA1D11_dge.tsv.gz
|       |   |   |-- timecourse_d7911_d11_MK_T21wtGATA1D11_dge.tsv.gz
|       |   |   `-- timecourse_d7911_d11_Myeloid_T21wtGATA1D11_dge.tsv.gz
|       |   `-- D7
|       |       |-- GSEA
|       |       |-- timecourse_d7911_d7_Ery_T21wtGATA1D7_dge.tsv.gz
|       |       |-- timecourse_d7911_d7_HPC_-_Ery_bias_T21wtGATA1D7_dge.tsv.gz
|       |       |-- timecourse_d7911_d7_HPC_-_MK_bias_1_T21wtGATA1D7_dge.tsv.gz
|       |       |-- timecourse_d7911_d7_HPC_-_MK_bias_2_T21wtGATA1D7_dge.tsv.gz
|       |       |-- timecourse_d7911_d7_HPCs_T21wtGATA1D7_dge.tsv.gz
|       |       |-- timecourse_d7911_d7_MK_T21wtGATA1D7_dge.tsv.gz
|       |       `-- timecourse_d7911_d7_Myeloid_T21wtGATA1D7_dge.tsv.gz
|       `-- timecourse-d7911-umap-renamed-cell-counts.tsv.gz
|-- seurat_objects
|   |-- data.info
|   |-- query_primary_plus_d7.RDS
|   |-- reference.RDS
|   |-- timecourse-d7-and-primary-renamed.RDS
|   |-- timecourse-d7-and-primary.RDS
|   |-- timecourse-d7911-renamed.RDS
|   `-- timecourse-d7911.RDS
`-- trajectory_analysis_notebooks
    |-- 07_prepare_d7911_seurat_for_anndata.R
    |-- 08_d7911_scv_create_anndata.ipynb
    |-- 09_d7911_scv_velocity_latent_time.ipynb
    |-- 10_d7911_cellrank_cell_fate.ipynb
    |-- 11_d7911_stalled_maturation_dge_analysis.R
    |-- 12_d7911_stalled_erythroid_maturation.ipynb
    |-- 12_d7911_stalled_megakaryrocyte_maturation.ipynb
    |-- 12_d7911_stalled_myeloid1_maturation.ipynb
    |-- 12_d7911_stalled_myeloid2_maturation.ipynb
    |-- 13_d7911_stalled_maturation_enrichment_analysis.Rmd
    `-- scVelo-latent-time-ridgeline-plots.Rmd
```


## Data availability 
- The 10x scRNA-Seq sequencing data generated in this study have been deposited in the NCBI Gene Expression Omnibus (GEO) under accession number GSE271399
- Complementary reference data 10x scRNA-Seq sequencing data was obtained from the [Calvanese et al., 2022](https://www.nature.com/articles/s41586-022-04571-x) study.

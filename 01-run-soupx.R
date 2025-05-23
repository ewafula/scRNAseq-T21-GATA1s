# Remove background RNA (contaminants) with SoupX

# Eric Wafula
# 10/20/2024

# Load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(SoupX))
suppressPackageStartupMessages(library(DoubletFinder))
suppressPackageStartupMessages(library(DropletUtils))

# function
run_soupx <- function(cellranger_dir, results_dir)
{
  sc <- SoupX::load10X(cellranger_dir)
  sc <- SoupX::autoEstCont(sc, soupQuantile = 0.75, tfidfMin = 0.5)
  results <- SoupX::adjustCounts(sc)
  DropletUtils::write10xCounts(results_dir, results, version='3', overwrite=TRUE)
}

# set up optparse options
option_list <- list(
  make_option(opt_str = "--project_dir", type = "character", default = NULL,
              help = "Project directory cellranger samples runs sub-directory",  
              metavar = "character")
)

# parse parameter options
opt <- parse_args(OptionParser(option_list = option_list))
project_dir <- opt$project_dir

# set seed for reproducibility
set.seed(123)

# establish base, data, and plots directories
root_dir <- here::here()
data_dir <- file.path(root_dir, "data")
plots_dir <- file.path(root_dir, "plots", "soupx")

# get samples names
samples <- list.files(here(data_dir, sample_type))

# remove contaminants - background (ambient) RNA
for (sample in samples) {
  # set input and output directories
  cellranger_dir <- file.path(data_dir, project_dir, sample, "outs")
  soupx_dir <- file.path(data_dir, project_dir, sample, "soupx")
  
  # create SoupX output directory if it doesn't exist
  if (!dir.exists(soupx_dir)) {
    dir.create(soupx_dir)
  }
  
  # remove sample contaminants with SoupX
  pdf(
    here::here(plots_dir, paste0(sample,"_soupx_contaminant.pdf")),
    height = 8,
    width = 11
  )
  run_soupx(cellranger_dir, soupx_dir)
  dev.off()
}


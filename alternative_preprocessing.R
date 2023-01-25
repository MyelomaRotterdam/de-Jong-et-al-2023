# --------------------------------------------------- #
# Analysis of scRNA sequencing data as described in de Jong et al. 2023
# Script by Madelon de Jong, optimized for speed by Mathijs Sanders
# Dept. of Hematology, Erasmus MC Cancer Institute, Rotterdam, the Netherlands
# --------------------------------------------------- #

# --------------------------------------------------- #
# READ-ME (IMPORTANT)
# In the paper, these three subsets are separated for individual analysis using selection of barcodes (available on Github), followed by reprocessing.
# Therefore, this function entails barcode selection to select out subsets of cells (neutrophils, monocytes, progenitors).
# Barcodes are available on the github.
# --------------------------------------------------- #

# --------------------------------------------------- #
# Loading libraries
# --------------------------------------------------- #

pacman::p_load('Seurat', 'future', 'parallel', 'cowplot', 'ggplot2', 'dplyr')

options(future.rng.onMisuse="ignore")
options(future.globals.maxSize = 1e4 * 1024^2)
options(future.fork.enable = TRUE)

# --------------------------------------------------- #
# Function for loading in objects
# --------------------------------------------------- #

loadSeuratObject <- function(idx, projectName, matLoc, bcLoc, selectionCriteria, details) {
  info <- structure(selectionCriteria[idx, ], names = colnames(selectionCriteria))
  currentDetails <- structure(details[idx, ], names = colnames(details))
  so <- Read10X(data.dir = matLoc[idx])
  so <- CreateSeuratObject(counts = so, min.cells = 3, min.features = 200, project = projectName)
  so <- AddMetaData(so, PercentageFeatureSet(so, pattern = "^MT-"), col.name = "percent.mito")
  so <- subset(x = so, subset = nFeature_RNA > info$nFeature_RNA_lower & nFeature_RNA < info$nFeature_RNA_higher & nCount_RNA > info$nCount_RNA_lower & nCount_RNA < info$nCount_RNA_higher & percent.mito < info$percent.mito)
  select.cells_so <- paste0(as.character(read.csv(file= bcLoc[idx], header=T, row.names = 1)[,1]))
  select.cells_so <- gsub("_.*", "", select.cells_so)
  so <- subset(x = so, cells = select.cells_so)
  so <- NormalizeData(object = so, normalization.method = "LogNormalize", scale.factor = 1e4)
  so <- FindVariableFeatures(object = so, selection.method = "vst", nfeatures = 2000)
  so$timepoint <- currentDetails$Timepoint
  so$patient <- currentDetails$Patient
  so$treatment <- currentDetails$Treatment
  so$conglom <- currentDetails$Conglomerate
  so$genetics <- currentDetails$Genetics
  so
}

# --------------------------------------------------- #
# Function for defining number of cores
# --------------------------------------------------- #

mcPredefined <- function(threads) {
  function(...) {
    mclapply(..., mc.cores = threads)
  }
}

# --------------------------------------------------- #
# Object information
# --------------------------------------------------- #

projectName <- 'Non-hematopoietic'
threads <- 10
matLoc <- readLines('~/samples.txt')
selectionCriteria <- read.table('~/selection.txt', header = TRUE)
details <- read.table('~/details.txt', header = TRUE)
bcLoc <- readLines('~/barcodes.txt')
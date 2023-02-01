# Script to estimate CD4+ and CD8+ T cells content in single-cell RNA-seq data from Gueguen P et al 2021.
# Gueguen P, Metoikidou C, Dupic T, Lawand M, Goudot C, Baulande S, Lameiras S, Lantz O, Girard N, Seguin-Givelet A, Lefevre M, Mora T, Walczak AM, Waterfall JJ, Amigorena S. Contribution of 
# resident and circulating precursors to tumor-infiltrating CD8+ T cell populations in lung cancer. Sci Immunol. 2021 Jan 29;6(55):eabd5778.

#Clean workspace
rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)

#---------- library -----------#
library(Seurat)
library(data.table)
library(dplyr)
#---------- library -----------#

# Load output of the cellranger pipeline from 10X, returning a unique molecular identified (UMI) count matrix for 9 patients.
setwd("/Functional_Heterogeneity/Public_Data/10X/Gueguen_2021_GSE162500/mappedData/SeuratInput")

folders <- list.files()

# Loop all over the folders and save the relevant data.
expressionList <- list()
for (folder in folders)
{
  expressionMatrix <- Read10X(data.dir = folder)
  colnames(expressionMatrix) <- paste0(folder, "_", colnames(expressionMatrix))
  
  expressionList[[folder]] <- expressionMatrix
}

# Merge the data together.
expressionMatrix <- do.call(cbind, expressionList)

############################ Preprocessing ###############################
# Initialize the Seurat object with the raw (non-normalized data).
seuratObject <- CreateSeuratObject(counts = expressionMatrix , min.cells = 3, min.features = 200, project = "CCIT")

# Create the data matrix.
dataMatrix <- GetAssayData(object = seuratObject)
features <- rownames(dataMatrix)
dataMatrix <- data.frame(as.matrix(dataMatrix))
dataMatrix <- cbind(Feature = features, dataMatrix)
dim(dataMatrix)

# Create the metadata matrix.
Cell <- colnames(dataMatrix)[-1]
groups <- strsplit(Cell, "_")
lengths <- sapply(groups, length)
for (i in 1:length(groups)) groups[[i]] <- groups[[i]][-lengths[[i]]]
groups <- sapply(groups, paste0, collapse = "_")

metadataMatrix <- data.frame(Cell = Cell, Patient = groups, stringsAsFactors = FALSE)

# Add new column for Histology Type, Sample Timing, Treatment and Dataset name.
metadataMatrix$Patient <- paste(metadataMatrix$Patient, "Gueguen_2021", sep ="_")
metadataMatrix$Histology <- "NSCLC"
metadataMatrix$HistologySubType <- "LUAD"
metadataMatrix$SampleTiming <- "Pre"
metadataMatrix$Treatment <- "Unknown"
metadataMatrix$CellType <- "NA"
metadataMatrix$Dataset <- "Gueguen_2021_GSE162500"

d_Gueguen_2021 <- dataMatrix
md_Gueguen_2021 <- metadataMatrix

# For merged data.
d_Gueguen_2021 <- d_Gueguen_2021[, -1]
dim(d_Gueguen_2021)

rownames(md_Gueguen_2021) <- md_Gueguen_2021[, 1]
md_Gueguen_2021 <- md_Gueguen_2021[, -1]

# Check.
identical(colnames(d_Gueguen_2021), rownames(md_Gueguen_2021))
rm(seuratObject)

d_s <- d_Gueguen_2021
d_ms <- md_Gueguen_2021

# Initialize the Seurat object with the raw (non-normalized data).
seuratObject <- CreateSeuratObject(counts= d_s, meta.data = d_ms, min.cells = 3, min.features = 200)

print(paste("Number of selected cells:", ncol(seuratObject@assays$RNA@counts)))
# "Number of selected cells: 58887"

print(paste("Number of selected genes:", nrow(seuratObject@assays$RNA@counts)))
# "Number of selected genes: 23928"

print(paste("Number of selected patients:", length(unique(seuratObject@meta.data$Patient)) ))
# "Number of selected patients: 9"

# Normalize using Seurat function SCTransform with batch_var=Patient to regress out latent variables.
allCell <- SCTransform(seuratObject, batch_var = "Patient", verbose = FALSE)

# Dimensionality reduction.
allCell <- RunPCA(object = allCell, npcs = 20, verbose = FALSE)
allCell <- RunTSNE(object = allCell, reduction = "pca",
                   dims = 1:20, check_duplicates = FALSE)

# Visualize dimensionality reduction results.
DimPlot(object = allCell, reduction = "tsne", group.by = "Patient")

# Get normalized data and metadata
d_s_norm <- as.data.frame(allCell@assays$SCT@data)
d_meta_norm <- as.data.frame(allCell@meta.data)

# Distribution of CD3D expression.
hist(as.numeric(d_s_norm['CD3D',]), breaks=10, main="Distribution of CD3D expression", xlab="Normalized counts")
abline(v=0, col="red", lwd=3)

# Assign the value of expression from the data matrix.
d_meta_norm$val_CD3D <- as.numeric(d_s_norm['CD3D',])

# Assign a logical value of according to the expression of the data matrix.
d_meta_norm$exp_CD3D <- as.numeric(d_s_norm['CD3D',])>= 0

# # Distribution of CD3E expression
hist(as.numeric(d_s_norm['CD3E',]), breaks=5, main="Distribution of CD3E expression", xlab="Normalized counts")
abline(v=0, col="red", lwd=3)

# Assign the value of expression from the data matrix.
d_meta_norm$val_CD3E <- as.numeric(d_s_norm['CD3E',])

# Assign a logical value of according to the expression of the data matrix.
d_meta_norm$exp_CD3E <- as.numeric(d_s_norm['CD3E',])>= 0

# # Distribution of CD3G expression.
hist(as.numeric(d_s_norm['CD3G',]), breaks=10, main="Distribution of CD3G expression", xlab="Normalized counts")
abline(v=0, col="red", lwd=3)

# Assign the value of expression from the data matrix.
d_meta_norm$val_CD3G <- as.numeric(d_s_norm['CD3G',])

# Assign a logical value of according to the expression of the data matrix.
d_meta_norm$exp_CD3G <- as.numeric(d_s_norm['CD3G',])>= 0

# Distribution of CD19 expression.
hist(as.numeric(d_s_norm['CD19',]), breaks=10, main="Distribution of CD19 expression", xlab="Normalized counts")
abline(v=0.6, col="red", lwd=3)

# Assign the value of expression from the data matrix.
d_meta_norm$val_CD19 <- as.numeric(d_s_norm['CD19',])

# Assign a logical value of according to the expression of the data matrix.
d_meta_norm$exp_CD19 <- as.numeric(d_s_norm['CD19',])> 0.6

# Distribution of MS4A1 expression.
hist(as.numeric(d_s_norm['MS4A1',]), breaks=10, main="Distribution of MS4A1 expression", xlab="Normalized counts")
abline(v=0.6, col="red", lwd=3)

# Assign the value of expression from the data matrix.
d_meta_norm$val_MS4A1 <- as.numeric(d_s_norm['MS4A1',])

# Assign a logical value of according to the expression of the data matrix.
d_meta_norm$exp_MS4A1 <- as.numeric(d_s_norm['MS4A1',])> 0.6

# Distribution of CD14 expression.
hist(as.numeric(d_s_norm['CD14',]), breaks=10, main="Distribution of CD14 expression", xlab="Normalized counts")
abline(v=0.6, col="red", lwd=3)

# Assign the value of expression from the data matrix.
d_meta_norm$val_CD14 <- as.numeric(d_s_norm['CD14',])

# Assign a logical value of according to the expression of the data matrix.
d_meta_norm$exp_CD14 <- as.numeric(d_s_norm['CD14',])> 0.6

# Distribution of CD4 expression.
hist(as.numeric(d_s_norm['CD4',]), breaks=10, main="Distribution of CD4 expression", xlab="Normalized counts")
abline(v=0.6, col="red", lwd=3)

# Assign the value of expression from the data matrix.
d_meta_norm$val_CD4 <- as.numeric(d_s_norm['CD4',])

# Assign a logical value of according to the expression of the data matrix.
d_meta_norm$exp_CD4 <- as.numeric(d_s_norm['CD4',])> 0.6

# Distribution of CD8A expression.
hist(as.numeric(d_s_norm['CD8A',]), breaks=10, main="Distribution of CD8A expression", xlab="Normalized counts")
abline(v=0.6, col="red", lwd=3)

# Assign the value of expression from the data matrix.
d_meta_norm$val_CD8A <- as.numeric(d_s_norm['CD8A',])

# Assign a logical value of according to the expression of the data matrix.
d_meta_norm$exp_CD8A <- as.numeric(d_s_norm['CD8A',])> 0.6

# Distribution of CD8B expression.
hist(as.numeric(d_s_norm['CD8B',]), breaks=10, main="Distribution of CD8B expression", xlab="Normalized counts")
abline(v=0.6, col="red", lwd=3)

# Assign the value of expression from the data matrix.
d_meta_norm$val_CD8B <- as.numeric(d_s_norm['CD8B',])

# Assign a logical value of according to the expression of the data matrix.
d_meta_norm$exp_CD8B <- as.numeric(d_s_norm['CD8B',])> 0.6

# Distribution of NCAM1 expression
hist(as.numeric(d_s_norm['NCAM1',]), breaks=10, main="Distribution of NCAM1 expression", xlab="Normalized counts")
abline(v=0.6, col="red", lwd=3)

# Assign the value of expression from the data matrix.
d_meta_norm$NCAM1 <- as.numeric(d_s_norm['NCAM1',])

# Assign a logical value of according to the expression of the data matrix.
d_meta_norm$exp_NCAM1 <- as.numeric(d_s_norm['NCAM1',])> 0.6

# Matrices only containing: CD3D+ OR CD3E+ OR CD3G+ cells.
index1 <- which((d_meta_norm$exp_CD3D == TRUE | d_meta_norm$exp_CD3E == TRUE | d_meta_norm$exp_CD3G == TRUE) & (d_meta_norm$exp_CD19 == FALSE & d_meta_norm$exp_MS4A1 == FALSE & d_meta_norm$exp_CD14 == FALSE))
length(index1)
# 56284

md_CD3D_Gueguen_2021_GSE162500 <- d_meta_norm[index1, ]
dim(md_CD3D_Gueguen_2021_GSE162500)
# 56284    32

d_CD3D_Gueguen_2021_GSE162500 <- subset(d_s_raw, select = rownames(md_CD3D_Gueguen_2021_GSE162500))
dim(d_CD3D_Gueguen_2021_GSE162500)
# 23928 56284

# Check
identical(colnames(d_CD3D_Gueguen_2021_GSE162500), rownames(md_CD3D_Gueguen_2021_GSE162500))

# Save raw data and metadata for CD3+ cells.
save(d_CD3D_Gueguen_2021_GSE162500, md_CD3D_Gueguen_2021_GSE162500, file="/Functional_Heterogeneity/Public_Data/10X/Gueguen_2021_GSE162500/processedData/d_md_CD3D_Gueguen_2021_GSE162500.Rdata")

# Matrices only containing: (CD3D+ OR CD3E+ OR CD3G+) AND CD4+ AND CD8A- AND CD8B- cells.
index2 <- which((d_meta_norm$exp_CD3D == TRUE | d_meta_norm$exp_CD3E == TRUE | d_meta_norm$exp_CD3G == TRUE) & d_meta_norm$exp_CD4 == TRUE & d_meta_norm$exp_CD8A == FALSE & d_meta_norm$exp_CD8B == FALSE & d_meta_norm$exp_CD19 == FALSE & d_meta_norm$exp_MS4A1 == FALSE & d_meta_norm$exp_CD14 == FALSE) 
length(index2)
# 5891

md_CD4_Gueguen_2021_GSE162500 <- d_meta_norm[index2, ]
dim(md_CD4_Gueguen_2021_GSE162500)
# 5891   32

d_CD4_Gueguen_2021_GSE162500 <- subset(d_s_raw, select = rownames(md_CD4_Gueguen_2021_GSE162500))
dim(d_CD4_Gueguen_2021_GSE162500)
# 23928  5891

# Check
identical(colnames(d_CD4_Gueguen_2021_GSE162500), rownames(md_CD4_Gueguen_2021_GSE162500))

# Save data and raw for CD4+ cells.
save(d_CD4_Gueguen_2021_GSE162500, md_CD4_Gueguen_2021_GSE162500, file="/Functional_Heterogeneity/Public_Data/10X/Gueguen_2021_GSE162500/processedData/d_md_CD4_Gueguen_2021_GSE162500.Rdata")

# Matrices only containing: (CD3D+ OR CD3E+ OR CD3G+) AND CD4- AND (CD8A+ OR CD8B+) cells.
index3 <- which((d_meta_norm$exp_CD3D == TRUE | d_meta_norm$exp_CD3E == TRUE | d_meta_norm$exp_CD3G == TRUE) & d_meta_norm$exp_CD4 == FALSE & d_meta_norm$exp_CD19 == FALSE & d_meta_norm$exp_MS4A1 == FALSE & d_meta_norm$exp_CD14 == FALSE & (d_meta_norm$exp_CD8A == TRUE | d_meta_norm$exp_CD8B == TRUE)) 
length(index3)
# 11679

md_CD8_Gueguen_2021_GSE162500 <- d_meta_norm[index3, ]
dim(md_CD8_Gueguen_2021_GSE162500)
# 11679    32

d_CD8_Gueguen_2021_GSE162500 <- subset(d_s_raw, select = rownames(md_CD8_Gueguen_2021_GSE162500))
dim(d_CD8_Gueguen_2021_GSE162500)
# 23928 11679

# Check
identical(colnames(d_CD8_Gueguen_2021_GSE162500), rownames(md_CD8_Gueguen_2021_GSE162500))

# Save data and raw for CD8+ cells.
save(d_CD8_Gueguen_2021_GSE162500, md_CD8_Gueguen_2021_GSE162500, file="/Functional_Heterogeneity/Public_Data/10X/Gueguen_2021_GSE162500/processedData/d_md_CD8_Gueguen_2021_GSE162500.Rdata")

# Script to estimate CD4+ and CD8+ T cells content in single-cell RNA-seq data from Bassez et al 2021.
# Bassez, A., Vos, H., Van Dyck, L. et al. A single-cell map of intratumoral changes during anti-PD1 treatment of patients with breast cancer. Nat Med 27, 820??832 (2021).

#Clean workspace
rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)

#---------- library -----------#
library(Seurat)
library(data.table)
library(dplyr)
#---------- library -----------#

# Read count data and metadata per individual patient were downloaded from http://biokey.lambrechtslab.org/

#---------- load data -----------#
cohortA <- readRDS("Functional heterogeneity/Public_Data/10X/Bassez_2021_EGAS00001004809/1863-counts_cells_cohort1.rds")

#---------- load metadata -----------#
md1 <- fread("E:/Analysis/Arianna/Functional heterogeneity/Public_Data/10X/Bassez_2021_EGAS00001004809/1872-BIOKEY_metaData_cohort1_web.csv", data.table = FALSE)

# Unique cell type from metadata
unique(md1$cellType)

# Select only T cells from metadata
idx1 <- which(md1$cellType == "T_cell")

# Subset metadata only for T cells
d_ms <- md1[idx1, ]

# Subset count data only for T cells
d_s <- cohortA[, d_ms$Cell]

d_s <- as.data.frame(as.matrix(d_s))

# Rename the metadata columns
colnames(d_ms)[colnames(d_ms) == 'patient_id'] <- 'Patient'
colnames(d_ms)[colnames(d_ms) == 'BC_type'] <- 'HistologySubType'
colnames(d_ms)[colnames(d_ms) == 'cellType'] <- 'CellType'

# Add new column for Histology Type, Sample Timing, Treatment and Dataset name
d_ms$Histology <- "BC"
d_ms$SampleTiming <- "Pre"
d_ms$Treatment <- "Unknown"
d_ms$Dataset <- "Bassez_2021_EGAS00001004809"

d_ms <- d_ms[, c("Cell", "Patient", "Histology", "HistologySubType", "SampleTiming", "Treatment", "CellType", "Dataset")]

#Preprocessing:  Basic cleaning:
rownames(d_ms) <- d_ms$Cell
d_ms <- d_ms[, -1]

# Check
identical(colnames(d_s), rownames(d_ms))

# Initialize the Seurat object with the raw (non-normalized data).
seuratObject  <- CreateSeuratObject(counts= d_s, meta.data = d_ms, min.cells = 3, min.features = 200)

print(paste("Number of selected cells:", ncol(seuratObject@assays$RNA@counts)))
# "Number of selected cells: 56968"

print(paste("Number of selected genes:", nrow(seuratObject@assays$RNA@counts)))
# "Number of selected genes: 20200"

print(paste("Number of selected patients:", length(unique(seuratObject@meta.data$Patient)) ))
# "Number of selected patients: 31"

# Normalize using Seurat function SCTransform with batch_var=Patient to regress out latent variables.
allCell <- SCTransform(seuratObject, batch_var = "Patient", verbose = FALSE)

# Get normalized data and metadata
d_s_norm <- as.data.frame(allCell@assays$SCT@data)
d_meta_norm <- as.data.frame(allCell@meta.data)

# Distribution of CD3D expression
hist(as.numeric(d_s_norm['CD3D',]), breaks=10, main="Distribution of CD3D expression", xlab="Normalized counts")
abline(v=0, col="red", lwd=3)

# Assign the value of expression from the data matrix
d_meta_norm$val_CD3D <- as.numeric(d_s_norm['CD3D',])

# Assign a logical value of according to the expression of the data matrix
d_meta_norm$exp_CD3D <- as.numeric(d_s_norm['CD3D',])>= 0

# # Distribution of CD3E expression
hist(as.numeric(d_s_norm['CD3E',]), breaks=5, main="Distribution of CD3E expression", xlab="Normalized counts")
abline(v=0, col="red", lwd=3)

# Assign the value of expression from the data matrix
d_meta_norm$val_CD3E <- as.numeric(d_s_norm['CD3E',])

# Assign a logical value of according to the expression of the data matrix
d_meta_norm$exp_CD3E <- as.numeric(d_s_norm['CD3E',])>= 0

# # Distribution of CD3G expression
hist(as.numeric(d_s_norm['CD3G',]), breaks=5, main="Distribution of CD3G expression", xlab="Normalized counts")
abline(v=0, col="red", lwd=3)

# Assign the value of expression from the data matrix
d_meta_norm$val_CD3G <- as.numeric(d_s_norm['CD3G',])

# Assign a logical value of according to the expression of the data matrix
d_meta_norm$exp_CD3G <- as.numeric(d_s_norm['CD3G',])>= 0

# Distribution of CD19 expression
hist(as.numeric(d_s_norm['CD19',]), breaks=10, main="Distribution of CD19 expression", xlab="Normalized counts")
abline(v=0.6, col="red", lwd=3)

# Assign the value of expression from the data matrix
d_meta_norm$val_CD19 <- as.numeric(d_s_norm['CD19',])

# Assign a logical value of according to the expression of the data matrix
d_meta_norm$exp_CD19 <- as.numeric(d_s_norm['CD19',])> 0.6

# Distribution of MS4A1 expression
hist(as.numeric(d_s_norm['MS4A1',]), breaks=10, main="Distribution of MS4A1 expression", xlab="Normalized counts")
abline(v=0.6, col="red", lwd=3)

# Assign the value of expression from the data matrix
d_meta_norm$val_MS4A1 <- as.numeric(d_s_norm['MS4A1',])

# Assign a logical value of according to the expression of the data matrix
d_meta_norm$exp_MS4A1 <- as.numeric(d_s_norm['MS4A1',])> 0.6

# Distribution of CD14 expression
hist(as.numeric(d_s_norm['CD14',]), breaks=10, main="Distribution of CD14 expression", xlab="Normalized counts")
abline(v=0.6, col="red", lwd=3)

# Assign the value of expression from the data matrix
d_meta_norm$val_CD14 <- as.numeric(d_s_norm['CD14',])

# Assign a logical value of according to the expression of the data matrix
d_meta_norm$exp_CD14 <- as.numeric(d_s_norm['CD14',])> 0.6

# Distribution of CD4 expression
hist(as.numeric(d_s_norm['CD4',]), breaks=10, main="Distribution of CD4 expression", xlab="Normalized counts")
abline(v=0.6, col="red", lwd=3)

# Assign the value of expression from the data matrix
d_meta_norm$val_CD4 <- as.numeric(d_s_norm['CD4',])

# Assign a logical value of according to the expression of the data matrix
d_meta_norm$exp_CD4 <- as.numeric(d_s_norm['CD4',])> 0.6

# Distribution of CD8A expression
hist(as.numeric(d_s_norm['CD8A',]), breaks=10, main="Distribution of CD8A expression", xlab="Normalized counts")
abline(v=0.6, col="red", lwd=3)

# Assign the value of expression from the data matrix
d_meta_norm$val_CD8A <- as.numeric(d_s_norm['CD8A',])

# Assign a logical value of according to the expression of the data matrix
d_meta_norm$exp_CD8A <- as.numeric(d_s_norm['CD8A',])> 0.6

# Distribution of CD8B expression
hist(as.numeric(d_s_norm['CD8B',]), breaks=10, main="Distribution of CD8B expression", xlab="Normalized counts")
abline(v=0.6, col="red", lwd=3)

# Assign the value of expression from the data matrix
d_meta_norm$val_CD8B <- as.numeric(d_s_norm['CD8B',])

# Assign a logical value of according to the expression of the data matrix
d_meta_norm$exp_CD8B <- as.numeric(d_s_norm['CD8B',])> 0.6

# Distribution of NCAM1 expression
hist(as.numeric(d_s_norm['NCAM1',]), breaks=10, main="Distribution of NCAM1 expression", xlab="Normalized counts")
abline(v=0.6, col="red", lwd=3)

# Assign the value of expression from the data matrix
d_meta_norm$NCAM1 <- as.numeric(d_s_norm['NCAM1',])

# Assign a logical value of according to the expression of the data matrix
d_meta_norm$exp_NCAM1 <- as.numeric(d_s_norm['NCAM1',])> 0.6

# Matrices only containing: CD3D+ OR CD3E+ OR CD3G+
index1 <- which((d_meta_norm$exp_CD3D == TRUE | d_meta_norm$exp_CD3E == TRUE | d_meta_norm$exp_CD3G == TRUE) & (d_meta_norm$exp_CD19 == FALSE & d_meta_norm$exp_MS4A1 == FALSE & d_meta_norm$exp_CD14 == FALSE))
length(index1)
#  54417

md_CD3D_Bassez_2021_EGAS00001004809 <- d_meta_norm[index1, ]
dim(md_CD3D_Bassez_2021_EGAS00001004809)
# 54417    32

d_CD3D_Bassez_2021_EGAS00001004809 <- subset(d_s_raw, select = rownames(md_CD3D_Bassez_2021_EGAS00001004809))
dim(d_CD3D_Bassez_2021_EGAS00001004809)
# 25288 54417

# Check
identical(colnames(d_CD3D_Bassez_2021_EGAS00001004809), rownames(md_CD3D_Bassez_2021_EGAS00001004809))

# Save raw data and metadata for CD3+ cells
save(d_CD3D_Bassez_2021_EGAS00001004809, md_CD3D_Bassez_2021_EGAS00001004809, file="/Functional_Heterogeneity/Public_Data/10X/Bassez_2021_EGAS00001004809/processedData/d_md_CD3D_Bassez_2021_EGAS00001004809.Rdata")

# Matrices only containing: (CD3D+ OR CD3E+ OR CD3G+) AND CD4+ AND CD8A- AND CD8B-
index2 <- which((d_meta_norm$exp_CD3D == TRUE | d_meta_norm$exp_CD3E == TRUE | d_meta_norm$exp_CD3G == TRUE) & d_meta_norm$exp_CD4 == TRUE & d_meta_norm$exp_CD8A == FALSE & d_meta_norm$exp_CD8B == FALSE & d_meta_norm$exp_CD19 == FALSE & d_meta_norm$exp_MS4A1 == FALSE & d_meta_norm$exp_CD14 == FALSE) 
length(index2)
# 10607

md_CD4_Bassez_2021_EGAS00001004809 <- d_meta_norm[index2, ]
dim(md_CD4_Bassez_2021_EGAS00001004809)
# 10607    32

d_CD4_Bassez_2021_EGAS00001004809 <- subset(d_s_raw, select = rownames(md_CD4_Bassez_2021_EGAS00001004809))
dim(d_CD4_Bassez_2021_EGAS00001004809)
# 25288 10607

# Check
identical(colnames(d_CD4_Bassez_2021_EGAS00001004809), rownames(md_CD4_Bassez_2021_EGAS00001004809))

# Save raw data and metadata for CD4+ cells
save(d_CD4_Bassez_2021_EGAS00001004809, md_CD4_Bassez_2021_EGAS00001004809, file="/Functional_Heterogeneity/Public_Data/10X/Bassez_2021_EGAS00001004809/processedData/d_md_CD4_Bassez_2021_EGAS00001004809.Rdata")

# Matrices only containing: (CD3D+ OR CD3E+ OR CD3G+) AND CD4- AND (CD8A+ OR CD8B+)
index3 <- which((d_meta_norm$exp_CD3D == TRUE | d_meta_norm$exp_CD3E == TRUE | d_meta_norm$exp_CD3G == TRUE) & d_meta_norm$exp_CD4 == FALSE & d_meta_norm$exp_CD19 == FALSE & d_meta_norm$exp_MS4A1 == FALSE & d_meta_norm$exp_CD14 == FALSE & (d_meta_norm$exp_CD8A == TRUE | d_meta_norm$exp_CD8B == TRUE)) 
length(index3)
# 21197

md_CD8_Bassez_2021_EGAS00001004809 <- d_meta_norm[index3, ]
dim(md_CD8_Bassez_2021_EGAS00001004809)
# 21197    32

d_CD8_Bassez_2021_EGAS00001004809 <- subset(d_s_raw, select = rownames(md_CD8_Bassez_2021_EGAS00001004809))
dim(d_CD8_Bassez_2021_EGAS00001004809)
# 25288 21197

# Check
identical(colnames(d_CD8_Bassez_2021_EGAS00001004809), rownames(md_CD8_Bassez_2021_EGAS00001004809))

# Save raw data and metadata for CD8+ cells
save(d_CD8_Bassez_2021_EGAS00001004809, md_CD8_Bassez_2021_EGAS00001004809, file="/Functional_Heterogeneity/Public_Data/10X/Bassez_2021_EGAS00001004809/processedData/d_md_CD8_Bassez_2021_EGAS00001004809.Rdata")


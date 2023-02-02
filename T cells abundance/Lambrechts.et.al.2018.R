# Script to estimate CD4+ and CD8+ T cells content in single-cell RNA-seq data from Lambrechts et al 2018.
# Lambrechts D, Wauters E, Boeckx B, Aibar S, Nittner D, Burton O, Bassez A, Decaluw�? H, Pircher A, Van den Eynde K, Weynand B, Verbeken E, De Leyn P, Liston A, Vansteenkiste J, Carmeliet P, 
# Aerts S, Thienpont B. Phenotype molding of stromal cells in the lung tumor microenvironment. Nat Med. 2018 Aug;24(8):1277-1289.

#Clean workspace
rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)

#---------- library -----------#
library(Seurat)
library(data.table)
library(dplyr)
#---------- library -----------#

# Read count data and metadata per individual patient were downloaded from http://biokey.lambrechtslab.org/.
setwd("/Functional_Heterogeneity/Public_Data/10X/Lambrechts_2018_E-MTAB-6149/mappedData/LC/export")

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
LC_MetadaData <- fread(file = "/Functional_Heterogeneity/Public_Data/10X/Lambrechts_2018_E-MTAB-6149/mappedData/LC/LC_metadata.csv.gz", data.table = FALSE)
LC_MetadaData$Cell <- paste0("LC_counts_", LC_MetadaData$Cell)

# Subset for Tumor cells.
LC_MetadaData <- filter(LC_MetadaData, CellFromTumor == TRUE)
dim(LC_MetadaData)

# Rename the columns.
colnames(LC_MetadaData)[colnames(LC_MetadaData) == 'PatientNumber'] <- 'Patient'
colnames(LC_MetadaData)[colnames(LC_MetadaData) == 'TumorType'] <- 'Histology'

LC_MetadaData$Patient <- paste0(LC_MetadaData$Histology, "_",  LC_MetadaData$Patient)
LC_MetadaData$Patient <- paste(LC_MetadaData$Patient, "Lambrechts_2018", sep ="_")

# Add new column for Histology Type, Sample Timing, Treatment and Dataset name
LC_MetadaData$Histology <- "NSCLC"
LC_MetadaData$HistologySubType <- "SCC"

index <- which(LC_MetadaData$Patient == "Lung_3_Lambrechts_2018")
LC_MetadaData$HistologySubType[index] <- "LUAD"

index <- which(LC_MetadaData$Patient == "Lung_4_Lambrechts_2018")
LC_MetadaData$HistologySubType[index] <- "LUAD"

index <- which(LC_MetadaData$Patient == "Lung_5_Lambrechts_2018")
LC_MetadaData$HistologySubType[index] <- "LC"

index <- which(LC_MetadaData$Patient == "Lung_6_Lambrechts_2018")
LC_MetadaData$HistologySubType[index] <- "LUAD"

index <- which(LC_MetadaData$Patient == "Lung_8_Lambrechts_2018")
LC_MetadaData$HistologySubType[index] <- "Pleio"

# Add new column for Histology Type, Sample Timing, Treatment and Dataset name.
LC_MetadaData$Histology <- "NSCLC"
LC_MetadaData$HistologySubType <- "SCC"

LC_MetadaData$SampleTiming <- "Pre"
LC_MetadaData$Treatment <- "Unknown"
LC_MetadaData$Dataset <- "Lambrechts_2018_E_MTAB_6149"
LC_MetadaData <- LC_MetadaData[, c("Cell", "Patient", "Histology", "HistologySubType", "SampleTiming", "Treatment", "CellType", "Dataset")]

rownames(LC_MetadaData) <- LC_MetadaData[ ,1]
LC_MetadaData <- LC_MetadaData[, -1]

dataMatrix <- select(dataMatrix, rownames(LC_MetadaData))

# Check.
identical(colnames(dataMatrix), rownames(LC_MetadaData))

# Subset for T cells.
idx <- which(LC_MetadaData$CellType == "T_cell")
LC_MetadaData <- LC_MetadaData[idx, ]
dataMatrix <- subset(dataMatrix, select = rownames(LC_MetadaData))
dim(dataMatrix)

# Check.
identical(colnames(dataMatrix), rownames(LC_MetadaData))

# Initialize the Seurat object with the raw (non-normalized data).
seuratObject <- CreateSeuratObject(counts= d_s, meta.data = d_ms, min.cells = 3, min.features = 200)

# Normalize using Seurat function SCTransform with batch_var=Patient to regress out latent variables.
allCell <- SCTransform(seuratObject, batch_var = "Patient", verbose = FALSE)

# Dimensionality reduction.
allCell <- RunPCA(object = allCell, npcs = 20, verbose = FALSE)
allCell <- RunTSNE(object = allCell, reduction = "pca",
                   dims = 1:20, check_duplicates = FALSE)

# Visualize dimensionality reduction results.
DimPlot(object = allCell, reduction = "tsne", group.by = "Patient")

# Get normalized data and metadata.
d_s_norm <- as.data.frame(allCell@assays$SCT@data)
d_meta_norm <- as.data.frame(allCell@meta.data)

# Distribution of CD3D expression.
hist(as.numeric(d_s_norm['CD3D',]), breaks=10, main="Distribution of CD3D expression", xlab="Normalized counts")
abline(v=0.6, col="red", lwd=3)

# Assign the value of expression from the data matrix.
d_meta_norm$val_CD3D <- as.numeric(d_s_norm['CD3D',])

# Assign a logical value of according to the expression of the data matrix.
d_meta_norm$exp_CD3D <- as.numeric(d_s_norm['CD3D',])> 0.6

# # Distribution of CD3E expression.
hist(as.numeric(d_s_norm['CD3E',]), breaks=5, main="Distribution of CD3E expression", xlab="Normalized counts")
abline(v=0.6, col="red", lwd=3)

# Assign the value of expression from the data matrix.
d_meta_norm$val_CD3E <- as.numeric(d_s_norm['CD3E',])

# Assign a logical value of according to the expression of the data matrix.
d_meta_norm$exp_CD3E <- as.numeric(d_s_norm['CD3E',])> 0.6

# # Distribution of CD3G expression.
hist(as.numeric(d_s_norm['CD3G',]), breaks=10, main="Distribution of CD3G expression", xlab="Normalized counts")
abline(v=0.6, col="red", lwd=3)

# Assign the value of expression from the data matrix.
d_meta_norm$val_CD3G <- as.numeric(d_s_norm['CD3G',])

# Assign a logical value of according to the expression of the data matrix.
d_meta_norm$exp_CD3G <- as.numeric(d_s_norm['CD3G',])> 0.6

# Distribution of CD14 expression
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

# Distribution of NCAM1 expression.
hist(as.numeric(d_s_norm['NCAM1',]), breaks=10, main="Distribution of NCAM1 expression", xlab="Normalized counts")
abline(v=0.6, col="red", lwd=3)

# Assign the value of expression from the data matrix.
d_meta_norm$NCAM1 <- as.numeric(d_s_norm['NCAM1',])

# Assign a logical value of according to the expression of the data matrix.
d_meta_norm$exp_NCAM1 <- as.numeric(d_s_norm['NCAM1',])> 0.6

# Matrices only containing: CD3D+ OR CD3E+ OR CD3G+ cells.
index1 <- which((d_meta_norm$exp_CD3D == TRUE | d_meta_norm$exp_CD3E == TRUE | d_meta_norm$exp_CD3G == TRUE) & (d_meta_norm$exp_CD19 == FALSE & d_meta_norm$exp_MS4A1 == FALSE & d_meta_norm$exp_CD14 == FALSE))
length(index1)

md_CD3D_Lambrechts_2018_E_MTAB_6149 <- d_meta_norm[index1, ]
dim(md_CD3D_Lambrechts_2018_E_MTAB_6149)

d_CD3D_Lambrechts_2018_E_MTAB_6149 <- subset(d_s_raw, select = rownames(md_CD3D_Lambrechts_2018_E_MTAB_6149))
dim(d_CD3D_Lambrechts_2018_E_MTAB_6149)

# Check.
identical(colnames(d_CD3D_Lambrechts_2018_E_MTAB_6149), rownames(md_CD3D_Lambrechts_2018_E_MTAB_6149))

# Save data and raw for CD3+ cells.
save(d_CD3D_Lambrechts_2018_E_MTAB_6149, md_CD3D_Lambrechts_2018_E_MTAB_6149, file="/Functional_Heterogeneity/Public_Data/10X/Lambrechts_2018_E-MTAB-6149/processedData/d_md_CD3D_Lambrechts_2018_E_MTAB_6149.Rdata")

# Matrices only containing: (CD3D+ OR CD3E+ OR CD3G+) AND CD4+ AND CD8A- AND CD8B- cells.
index2 <- which((d_meta_norm$exp_CD3D == TRUE | d_meta_norm$exp_CD3E == TRUE | d_meta_norm$exp_CD3G == TRUE) & d_meta_norm$exp_CD4 == TRUE & d_meta_norm$exp_CD8A == FALSE & d_meta_norm$exp_CD8B == FALSE & d_meta_norm$exp_CD19 == FALSE & d_meta_norm$exp_MS4A1 == FALSE & d_meta_norm$exp_CD14 == FALSE) 
length(index2)

md_CD4_Lambrechts_2018_E_MTAB_6149 <- d_meta_norm[index2, ]
dim(md_CD4_Lambrechts_2018_E_MTAB_6149)

d_CD4_Lambrechts_2018_E_MTAB_6149 <- subset(d_s_raw, select = rownames(md_CD4_Lambrechts_2018_E_MTAB_6149))
dim(d_CD4_Lambrechts_2018_E_MTAB_6149)

# Check.
identical(colnames(d_CD4_Lambrechts_2018_E_MTAB_6149), rownames(md_CD4_Lambrechts_2018_E_MTAB_6149))

# Save data and raw for CD4+ cells.
save(d_CD4_Lambrechts_2018_E_MTAB_6149, md_CD4_Lambrechts_2018_E_MTAB_6149, file="/Functional_Heterogeneity/Public_Data/10X/Lambrechts_2018_E-MTAB-6149/processedData/d_md_CD4_Lambrechts_2018_E_MTAB_6149.Rdata")

# Matrices only containing: (CD3D+ OR CD3E+ OR CD3G+) AND CD4- AND (CD8A+ OR CD8B+) cells.
index3 <- which((d_meta_norm$exp_CD3D == TRUE | d_meta_norm$exp_CD3E == TRUE | d_meta_norm$exp_CD3G == TRUE) & d_meta_norm$exp_CD4 == FALSE & d_meta_norm$exp_CD19 == FALSE & d_meta_norm$exp_MS4A1 == FALSE & d_meta_norm$exp_CD14 == FALSE & (d_meta_norm$exp_CD8A == TRUE | d_meta_norm$exp_CD8B == TRUE)) 
length(index3)

md_CD8_Lambrechts_2018_E_MTAB_6149 <- d_meta_norm[index3, ]
dim(md_CD8_Lambrechts_2018_E_MTAB_6149)

d_CD8_Lambrechts_2018_E_MTAB_6149 <- subset(d_s_raw, select = rownames(md_CD8_Lambrechts_2018_E_MTAB_6149))
dim(d_CD8_Lambrechts_2018_E_MTAB_6149)

# Check..
identical(colnames(d_CD8_Lambrechts_2018_E_MTAB_6149), rownames(md_CD8_Lambrechts_2018_E_MTAB_6149))

# Save data and raw for CD8+ cells
save(d_CD8_Lambrechts_2018_E_MTAB_6149, md_CD8_Lambrechts_2018_E_MTAB_6149, file="/Functional_Heterogeneity/Public_Data/10X/Lambrechts_2018_E-MTAB-6149/processedData/d_md_CD8_Lambrechts_2018_E_MTAB_6149.Rdata")

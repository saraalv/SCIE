## Preparación de los datos 

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DT)
library(GEOquery)

##------------------------------------------------------------
# 1. Obtención de los datos -- TCGA
##------------------------------------------------------------

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
)
GDCdownload(query = query)
data <- GDCprepare(query = query)

save(data, file = "C:/Users/Suruxx/Documents/TFM/raw data/tcga_rawdata1.RData")

##-------------------------------------------------------------
# 2. Pre-procesado de los datos -- TCGA
##-------------------------------------------------------------

## obtenemos un objeto tipo RangedSummarizedExperiment 
## guardamos los datos de expresión en un nuevo objeto y hacemos únicos los 
## nombres de los genes 
tcga_eset <- assay(data)

rownames(tcga_eset) <-  rowData(data)$gene_name
rownames(tcga_eset) <- make.unique(rownames(tcga_eset))

## creamos un nuevo objeto que contendrá la metadata

tcga_metadata <- subset(colData(data), select = c("barcode", "sample_type", "classification_of_tumor", "race", "gender", "ethnicity",
                                                  "primary_diagnosis"))

## guardamos los objetos en un archivo RData
save(tcga_eset, tcga_metadata, file = "C:/Users/Suruxx/Documents/TFM/raw data/tcga_data.RData")


##-----------------------------------------------------------------
# 3. Obtención de los datos -- AURORA 
##-----------------------------------------------------------------

# Cargamos el archivo descargado de GEO bajo el identificador GSE209998
# leemos los datos de expresión y hacemos únicos los nombres de los genes 
expr <- read.table("C:/Users/Suruxx/Documents/TFM/raw data/GSE209998_AUR_129_raw_counts.txt", sep = "\t", header = TRUE, quote = "")
expr$X <- make.unique(expr$X)
rownames(expr) <- expr$X
expr <- expr[, -1]

# Leemos los datos clínicos
clin <- read.csv2("C:/Users/Suruxx/Documents/TFM/raw data/clindata_aurora.csv")

##---------------------------------------------------------------
# 4. Pre-procesado de datos -- AURORA
##---------------------------------------------------------------

# Ajustamos los nombres de los IDs
names(expr) <- sub("^(.*?)\\.R\\..*", "\\1", names(expr))

# Filtramos las muestras en los datos clínicos 
expr <- expr[, names(expr) %in% clin$BCR.Portion.barcode]
clin <- clin[clin$BCR.Portion.barcode %in% names(expr), ]

# Alineamos los datos clínicos con los datos de expresión
colnames(expr) <- trimws(colnames(expr))
clin <- clin[match(colnames(expr), clin$BCR.Portion.barcode), ]

# Guardamos los datos 
save(expr, clin, file = "C:/Users/Suruxx/Documents/TFM/raw data/aurora_data.RData" )
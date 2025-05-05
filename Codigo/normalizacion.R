## Normalización 

library(edgeR)
library(WGCNA)
library(ggplot2)


#cargamos los datos 
load(file = "C:/Users/Suruxx/Documents/TFM/raw data/tcga_data.RData")
load(file = "C:/Users/Suruxx/Documents/TFM/raw data/aurora_data.RData")


##-------------------------------------------------------------
# 1. Detección de muestras y genes outliers 
##-------------------------------------------------------------

#detectamos genes atípicos 

gsg <- goodSamplesGenes(t(tcga_eset)) #tranponemos la matriz ya que la función pide que los genes sean las columnas y las filas las muestras
summary(gsg)
gsg$allOK #si sale TRUE significa que todo está bien, sino hay outliers

table(gsg$goodGenes)    # todos los genes pasan el control? --> los que tienen FALSE tienen demasiados missing values
table(gsg$goodSamples) #todas las muestras pasan el control?

gsg_aur <- goodSamplesGenes(t(expr))
gsg_aur$allOK
table(gsg_aur$goodGenes)
table(gsg_aur$goodSamples)

# eliminamos los genes detectados como outliers
tcga_eset <- tcga_eset[gsg$goodGenes == TRUE,]
aurora_eset <- expr[gsg_aur$goodGenes == TRUE,]

# detectamos muestras atípicas -  clustering jerarquico - método 1

htree <- hclust(dist(t(tcga_eset)), method = "average")
plot(htree)

ind <- which(htree$height > 5000000) # en el dendograma vemos que las muestras con altura > 5e6 son outliers
outl <- htree$labels[ind]

htree_aur <- hclust(dist(t(aurora_eset)), method = "average")
plot(htree_aur)
outl_aur <- c("AUR.AD9H.TTM1.A.2.1")

# Análisis de componentes principales - método 2 

pca <- prcomp(t(tcga_eset))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

pca_aur <- prcomp(t(aurora_eset))
pca.dat <- pca_aur$x

pca.var <- pca_aur$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

#eliminamos muestras detectadas como ouliers 
outl <- c(outl, "TCGA-EW-A1J6-01A-11R-A13Q-07")
outl_aur <- c(outl_aur, "AUR.AG12.TTM1.A.1.1","AUR.AER8.TTM2.A.1.1", "AUR.AFUK.TTM1.A.1.0")
tcga.subset <- tcga_eset[,!(colnames(tcga_eset) %in% outl)]
aurora.subset <- aurora_eset[,!colnames(aurora_eset) %in% outl_aur]

##-------------------------------------------------------------------
# 2. Eliminación de genes poco expresados
##-------------------------------------------------------------------

# TCGA
counts.CPM <- cpm(tcga.subset)
thresh <- counts.CPM >0.5
keep <- rowSums(thresh) >=2
counts.keep.tcga <- tcga.subset[keep,]

# AURORA
counts.CPM <- cpm(aurora.subset)
thresh <- counts.CPM >0.5
keep <- rowSums(thresh) >=2
counts.keep.aur <- aurora.subset[keep,]



## -----------------------------------------------------------
# 3. Normalización
## -----------------------------------------------------------

dge <- DGEList(counts.keep.tcga)       # guardamos los datos en un objeto dge
dge_norm <- calcNormFactors(dge)  # recalculamos los factores de normalización 
tcga_data <- cpm(dge_norm)

dge.aur <- DGEList(counts.keep.aur)
dge_norm <- calcNormFactors(dge.aur)
aurora_data <- cpm(dge_norm)

## --------------------------------------------------------------
# 4. Filtraje de las muestras en los metadatos
## --------------------------------------------------------------

aurora_data <- aurora_data[, colnames(aurora_data) %in% clin$BCR.Portion.barcode]
clin <- clin[clin$BCR.Portion.barcode %in% colnames(aurora_data),]

tcga_data <- tcga_data[, colnames(tcga_data) %in% tcga_metadata$barcode]
tcga_metadata <- tcga_metadata[tcga_metadata$barcode %in% colnames(tcga_data), ]


# guardamos los datos 
aurora_metadata <- clin
save(tcga_data, tcga_metadata, file = "C:/Users/Suruxx/Documents/TFM/data/tcga_dataset.RData")
save(aurora_data, aurora_metadata, file = "C:/Users/Suruxx/Documents/TFM/data/aurora_dataset.RData" )

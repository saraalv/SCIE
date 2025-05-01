#-----------------------------------------------------------------------------
## CREACIÓN DE LA RED DE COEXPRESIÓN -- AURORA  
#-----------------------------------------------------------------------------

# Cargamos los datos y las librerías necesarias 

library(WGCNA)
library(ggplot2)
library(gridExtra)
library(GSVA)
library(pheatmap)

load(file = "C:/Users/Suruxx/Documents/TFM/data/aurora_dataset.RData")


## -------------------------------------------------------------------------
# 1. Elección del umbral 
## -------------------------------------------------------------------------

power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# utilizamos la funcion pickSoftThreshold del paquete WGCNA
# que analiza la similitud con una red topologica 

sft <- pickSoftThreshold(t(aurora_data),
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices

## Visualizamos los resultados para elegir umbral 

## nos interesa un R^2 alto 

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.90, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

## y una conectividad media baja

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)

## -------------------------------------------------------------------------
# 2. Creación de la red 
## -------------------------------------------------------------------------


umbral <- 6
temp_cor <- cor
cor <- WGCNA::cor

bwnet <- blockwiseModules(t(aurora_data),
                          maxBlockSize = 10000,
                          TOMType = "signed",
                          power = umbral,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)
cor <- temp_cor

save(bwnet, file = "C:/Users/Suruxx/Documents/TFM/data/AURORA_red.RData")


## obtenemos el eigengene de cada módulo 
module_eigengenes <- bwnet$MEs 


## y el numero de genes que contiene cada módulo 
table(bwnet$colors)


## Representamos el dendograma y los módulos antes y depués del merging

plotDendroAndColors(bwnet$dendrograms[[1]], 
                    cbind(bwnet$unmergedColors[bwnet$blockGenes[[1]]], bwnet$colors[bwnet$blockGenes[[1]]]),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05,
                    main = "Cluster Dendogram block 1")

plotDendroAndColors(bwnet$dendrograms[[2]], 
                    cbind(bwnet$unmergedColors[bwnet$blockGenes[[2]]], bwnet$colors[bwnet$blockGenes[[2]]]),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05,
                    main = "Cluster Dendogram block 2")

plotDendroAndColors(bwnet$dendrograms[[3]], 
                    cbind(bwnet$unmergedColors[bwnet$blockGenes[[3]]], bwnet$colors[bwnet$blockGenes[[3]]]),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05,
                    main = "Cluster Dendogram block 3")
plotDendroAndColors(bwnet$dendrograms[[4]], 
                    cbind(bwnet$unmergedColors[bwnet$blockGenes[[4]]], bwnet$colors[bwnet$blockGenes[[4]]]),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05,
                    main = "Cluster Dendogram block 4")


## ------------------------------------------------------------------------
# 3. Obtención de los módulos de interes
## ------------------------------------------------------------------------

## creamos un dataframe con cada uno de los genes y el modulo al que pertenece 
modules_df <- as.data.frame(bwnet$colors)
modules_df <- cbind(names(bwnet$colors),modules_df)
colnames(modules_df) <- c("gene_id", "module")

## Mantenemos los modulos que contengan al menos 20 genes de firmas stem o inmune

## cargamos la información de las firmas
load("C:/Users/Suruxx/Documents/TFM/data/firmasgenicas.RData")

modules <- names(table(modules_df$module))
modulos_interes20 <- c()

inmune <- c(firmas[[1]],firmas[[2]], firmas[[3]], firmas[[4]])
stem <- c(firmas[[5]],firmas[[6]], firmas[[7]], firmas[[8]])

for (mod in modules) {
  f <- stem %in% subset(modules_df$gene_id, modules_df$module == mod)
  if (TRUE %in% f) {
    if(table(f)[[2]] >= 20){
      if( mod %in% modulos_interes20){
        
      }else{
        modulos_interes20 <- c(modulos_interes20, mod)
      }
    }
    
  }
  
  f <- inmune %in% subset(modules_df$gene_id, modules_df$module == mod)
  if (TRUE %in% f){
    if(table(f)[[2]] >= 20){
      if (mod %in% modulos_interes20){
        
      }else{
        modulos_interes20 <- c(modulos_interes20, mod)
      }
    }
  }
  
}


## -------------------------------------------------------------------------
# 4. Correlación entre los modúlos y las diferentes firmas 
## -------------------------------------------------------------------------

## calculamos la puntuacion gsva 

gsvapar <- gsvaParam(aurora_data, firmas)
gsva_aurora <- gsva(gsvapar, verbose = FALSE)

## obtenemos los datos de expresión por cada módulo
expr_bisque4 <- aurora_data[subset(modules_df$gene_id, modules_df$module == "bisque4"),]
expr_black <- aurora_data[subset(modules_df$gene_id, modules_df$module == "black"),]
expr_blue <- aurora_data[subset(modules_df$gene_id, modules_df$module == "blue"),]
expr_blue4 <- aurora_data[subset(modules_df$gene_id, modules_df$module == "blue4"),]
expr_brown <- aurora_data[subset(modules_df$gene_id, modules_df$module == "brown"),]
expr_darkgreen <- aurora_data[subset(modules_df$gene_id, modules_df$module == "darkgreen"),]
expr_darkmagenta <- aurora_data[subset(modules_df$gene_id, modules_df$module == "darkmagenta"),]
expr_darkred <- aurora_data[subset(modules_df$gene_id, modules_df$module == "darkred"),]
expr_darkseagreen4 <- aurora_data[subset(modules_df$gene_id, modules_df$module == "darkseagreen4"),]
expr_darkturquoise <- aurora_data[subset(modules_df$gene_id, modules_df$module == "darkturquoise"),]
expr_green <- aurora_data[subset(modules_df$gene_id, modules_df$module == "green"),]
expr_greenyellow <- aurora_data[subset(modules_df$gene_id, modules_df$module == "greenyellow"),]
expr_magenta <- aurora_data[subset(modules_df$gene_id, modules_df$module == "magenta"),]
expr_maroon <- aurora_data[subset(modules_df$gene_id, modules_df$module == "maroon"),]
expr_midnightblue <- aurora_data[subset(modules_df$gene_id, modules_df$module == "midnightblue"),]
expr_pink <- aurora_data[subset(modules_df$gene_id, modules_df$module == "pink"),]
expr_tan <- aurora_data[subset(modules_df$gene_id, modules_df$module == "tan"),]
expr_turquoise <- aurora_data[subset(modules_df$gene_id, modules_df$module == "turquoise"),]
expr_yellow <- aurora_data[subset(modules_df$gene_id, modules_df$module == "yellow"),]
expr_yellowgreen <- aurora_data[subset(modules_df$gene_id, modules_df$module == "yellowgreen"),]



# se repiten los siguientes pasos para cada módulo: 

## calculamos la media de expresión por muestra 
score_yellowgreen <- colMeans(expr_yellowgreen)

## Nos aseguramos que las muestras coincidan 
common_samples <- intersect(names(score_yellowgreen), colnames(gsva_tcga))

## calculamos la correlación con cada firma 
cor_isds_yellowgreen <- cor.test(score_yellowgreen[common_samples], gsva_tcga["ISDS", common_samples])
cor_mp17_yellowgreen <- cor.test(score_yellowgreen[common_samples], gsva_tcga["MP17", common_samples])
cor_mp18_yellowgreen <- cor.test(score_yellowgreen[common_samples], gsva_tcga["MP18", common_samples])
cor_kegg_yellowgreen <- cor.test(score_yellowgreen[common_samples], gsva_tcga["KEGG", common_samples])

cor_assou_yellowgreen <- cor.test(score_yellowgreen[common_samples], gsva_tcga["assou", common_samples])
cor_wong_yellowgreen <- cor.test(score_yellowgreen[common_samples], gsva_tcga["wong", common_samples])
cor_plurinet_yellowgreen <- cor.test(score_yellowgreen[common_samples], gsva_tcga["plurinet", common_samples])
cor_benporath_yellowgreen <- cor.test(score_yellowgreen[common_samples], gsva_tcga["benporath", common_samples])


## Matriz y heatmap 

ord_yellowgreen <- order(score_yellowgreen, decreasing = TRUE)
mat_yellowgreen <- rbind(
  ISDS = gsva_tcga["ISDS", ord_yellowgreen],
  MP17  = gsva_tcga["MP17", ord_yellowgreen],
  MP18  = gsva_tcga["MP18", ord_yellowgreen],
  KEGG  = gsva_tcga["KEGG", ord_yellowgreen],
  assou  = gsva_tcga["assou", ord_yellowgreen],
  wong  = gsva_tcga["wong", ord_yellowgreen],
  plurinet  = gsva_tcga["plurinet", ord_yellowgreen],
  benporath  = gsva_tcga["benporath", ord_yellowgreen]
)
mat_yellowgreen_scaled <- t(scale(t(mat_yellowgreen)))

labels_yellowgreen <- c(
  sprintf("ISDS signature\nr = %.2f, p = %.1e", cor_isds_yellowgreen$estimate, cor_isds_yellowgreen$p.value),
  sprintf("MP17 signature\nr = %.2f, p = %.1e", cor_mp17_yellowgreen$estimate, cor_mp17_yellowgreen$p.value),
  sprintf("MP18 signature\nr = %.2f, p = %.1e", cor_mp18_yellowgreen$estimate, cor_mp18_yellowgreen$p.value),
  sprintf("KEGG signature\nr = %.2f, p = %.1e", cor_kegg_yellowgreen$estimate, cor_kegg_yellowgreen$p.value),
  sprintf("Assou signature\nr = %.2f, p = %.1e", cor_assou_yellowgreen$estimate, cor_assou_yellowgreen$p.value),
  sprintf("Wong signature\nr = %.2f, p = %.1e", cor_wong_yellowgreen$estimate, cor_wong_yellowgreen$p.value),
  sprintf("Plurinet signature\nr = %.2f, p = %.1e", cor_plurinet_yellowgreen$estimate, cor_plurinet_yellowgreen$p.value),
  sprintf("Benporath signature\nr = %.2f, p = %.1e", cor_benporath_yellowgreen$estimate, cor_benporath_yellowgreen$p.value)
)

# ----------------------------
# 9. Colores y plots
# ----------------------------

heat_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# 📊 YELLOWGREEN
pheatmap(mat_yellowgreen_scaled,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = FALSE,
         show_rownames = TRUE,
         labels_row = labels_yellowgreen,
         fontsize_row = 9,
         color = heat_colors,
         main = "YELLOWGREEN MODULE:\n CORRELATION WITH GENETIC SIGNATURES")


## obtenemos los module eigengenes de los modulos de interes 

modulos_interesCOR <- c("blue4", "brown", "darkmagenta","darkseagreen4","greenyellow","maroon","yellow","yellowgreen")
MEmodules <- paste0("ME",modulos_interesCOR)

eigen_subset <- bwnet$MEs[,MEmodules]


## Correlacion entre los modulos de interés

##dendograma modulos de interés
plotEigengeneNetworks(eigen_subset,
                      names(eigen_subset),
                      excludeGrey = TRUE, greyLabel = "grey",
                      plotDendrograms = TRUE, plotHeatmaps = FALSE,
                      colorLabels = TRUE, signed = TRUE,
                      
                      plotAdjacency = FALSE)  # representar correlacion en vez de adyacencia 

## heatmap modulos de interés
plotEigengeneNetworks(eigen_subset,
                      names(eigen_subset),
                      excludeGrey = TRUE, greyLabel = "grey",
                      plotDendrograms = FALSE, plotHeatmaps = TRUE,
                      colorLabels = TRUE, signed = TRUE,
                      plotAdjacency = TRUE)  # representar correlacion en vez de adyacencia 

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
# 1. SEPARACIÓN DEL DATASET
## -------------------------------------------------------------------------



AURORA_ids <- subset(aurora_metadata$BCR.Portion.barcode, aurora_metadata$Sample.Type == "Primary")
aurora_primarios <- aurora_data[,AURORA_ids]
aurora_metastasis <- aurora_data[, !(colnames(aurora_data) %in% AURORA_ids)]

## --------------------------------------------------------------------------
# 2. ELECCIÓN DEL UMBRAL  Elección de umbral 
## --------------------------------------------------------------------------

power <- c(1:30)

sft <- pickSoftThreshold(t(aurora_metastasis),
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft.data <- sft$fitIndices

## Visualizamos los resultados para elegir umbral 

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.85, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

## --------------------------------------------------------------------
# 3. CREACIÓN DE LA RED 
## --------------------------------------------------------------------

umbral <- 14
temp_cor <- cor
cor <- WGCNA::cor

bwnet <- blockwiseModules(t(aurora_metastasis),
                          maxBlockSize = 10000,
                          TOMType = "signed",
                          power = umbral,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)
cor <- temp_cor

save(bwnet, file = "C:/Users/Suruxx/Documents/TFM/data/AURORAmet_red14.RData")

plotDendroAndColors(bwnet$dendrograms[[1]], 
                    cbind(bwnet$unmergedColors[bwnet$blockGenes[[1]]], bwnet$colors[bwnet$blockGenes[[1]]]),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05,
                    main = "Cluster Dendogram block 1 AURORA metastasis")

plotDendroAndColors(bwnet$dendrograms[[2]], 
                    cbind(bwnet$unmergedColors[bwnet$blockGenes[[2]]], bwnet$colors[bwnet$blockGenes[[2]]]),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05,
                    main = "Cluster Dendogram block 2 AURORA metastasis")

plotDendroAndColors(bwnet$dendrograms[[3]], 
                    cbind(bwnet$unmergedColors[bwnet$blockGenes[[3]]], bwnet$colors[bwnet$blockGenes[[3]]]),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05,
                    main = "Cluster Dendogram block 3 AURORA metastasis")
plotDendroAndColors(bwnet$dendrograms[[4]], 
                    cbind(bwnet$unmergedColors[bwnet$blockGenes[[4]]], bwnet$colors[bwnet$blockGenes[[4]]]),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05,
                    main = "Cluster Dendogram block 4 AURORA metastasis")

## -------------------------------------------------------------------------
# 4 . OBTENCIÓN DE LOS MÓDULOS DE INTERÉS 
## -------------------------------------------------------------------------

modules_df_met <- as.data.frame(bwnet$colors)
modules_df_met <- cbind(names(bwnet$colors),modules_df_met)
colnames(modules_df_met) <- c("gene_id", "module")

save(modules_df_met, file = "C:/Users/Suruxx/Documents/TFM/data/df_modulos_aurora_met.RData")

load("C:/Users/Suruxx/Documents/TFM/data/firmasgenicas.RData")

modules <- names(table(modules_df_met$module))
modulos_interes_met <- c()

inmune <- c(firmas[[1]],firmas[[2]], firmas[[3]], firmas[[4]])
stem <- c(firmas[[5]],firmas[[6]], firmas[[7]], firmas[[8]])

for (mod in modules) {
  f <- stem %in% subset(modules_df_met$gene_id, modules_df_met$module == mod)
  if (TRUE %in% f) {
    if(table(f)[[2]] >= 20){
      if( mod %in% modulos_interes_met){
        
      }else{
        modulos_interes_met <- c(modulos_interes_met, mod)
      }
    }
    
  }
  
  f <- inmune %in% subset(modules_df_met$gene_id, modules_df_met$module == mod)
  if (TRUE %in% f){
    if(table(f)[[2]] >= 20){
      if (mod %in% modulos_interes_met){
        
      }else{
        modulos_interes_met <- c(modulos_interes_met, mod)
      }
    }
  }
  
}

## --------------------------------------------------------------------
# 5. CORRELACIÓN ENTRE MÓDULOS Y FIRMAS 
## --------------------------------------------------------------------

## calculamos la puntuacion gsva 

gsvapar <- gsvaParam(aurora_metastasis, firmas)
gsva_aurora_met <- gsva(gsvapar, verbose = FALSE)

## obtenemos los datos de expresión por cada módulo

expr_blue_m <- aurora_metastasis[subset(modules_df_met$gene_id, modules_df_met$module == "blue"),]
expr_cyan_m <- aurora_metastasis[subset(modules_df_met$gene_id, modules_df_met$module == "cyan"),]
expr_midnightblue_m <- aurora_metastasis[subset(modules_df_met$gene_id, modules_df_met$module == "midnightblue"),]
expr_pink_m <- aurora_metastasis[subset(modules_df_met$gene_id, modules_df_met$module == "pink"),]
expr_red_m <- aurora_metastasis[subset(modules_df_met$gene_id, modules_df_met$module == "red"),]
expr_salmon_m <- aurora_metastasis[subset(modules_df_met$gene_id, modules_df_met$module == "salmon"),]
expr_turquoise_m <- aurora_metastasis[subset(modules_df_met$gene_id, modules_df_met$module == "turquoise"),]
expr_yellow_m <- aurora_metastasis[subset(modules_df_met$gene_id, modules_df_met$module == "yellow"),]



# se repiten los siguientes pasos para cada módulo: 

## calculamos la media de expresión por muestra 
score_turquoise <- colMeans(expr_turquoise_m)

## Nos aseguramos que las muestras coincidan 
common_samples <- intersect(names(score_turquoise), colnames(gsva_aurora_met))

## calculamos la correlación con cada firma 
cor_isds_turquoise <- cor.test(score_turquoise[common_samples], gsva_aurora_met["ISDS", common_samples])
cor_mp17_turquoise <- cor.test(score_turquoise[common_samples], gsva_aurora_met["MP17", common_samples])
cor_mp18_turquoise <- cor.test(score_turquoise[common_samples], gsva_aurora_met["MP18", common_samples])
cor_kegg_turquoise <- cor.test(score_turquoise[common_samples], gsva_aurora_met["KEGG", common_samples])

cor_assou_turquoise <- cor.test(score_turquoise[common_samples], gsva_aurora_met["assou", common_samples])
cor_wong_turquoise <- cor.test(score_turquoise[common_samples], gsva_aurora_met["wong", common_samples])
cor_plurinet_turquoise <- cor.test(score_turquoise[common_samples], gsva_aurora_met["plurinet", common_samples])
cor_benporath_turquoise <- cor.test(score_turquoise[common_samples], gsva_aurora_met["benporath", common_samples])


## Matriz y heatmap 

ord_turquoise <- order(score_turquoise, decreasing = TRUE)
mat_turquoise <- rbind(
  ISDS = gsva_aurora_met["ISDS", ord_turquoise],
  MP17  = gsva_aurora_met["MP17", ord_turquoise],
  MP18  = gsva_aurora_met["MP18", ord_turquoise],
  KEGG  = gsva_aurora_met["KEGG", ord_turquoise],
  assou  = gsva_aurora_met["assou", ord_turquoise],
  wong  = gsva_aurora_met["wong", ord_turquoise],
  plurinet  = gsva_aurora_met["plurinet", ord_turquoise],
  benporath  = gsva_aurora_met["benporath", ord_turquoise]
)
mat_turquoise_scaled <- t(scale(t(mat_turquoise)))

labels_turquoise <- c(
  sprintf("ISDS signature\nr = %.2f, p = %.1e", cor_isds_turquoise$estimate, cor_isds_turquoise$p.value),
  sprintf("MP17 signature\nr = %.2f, p = %.1e", cor_mp17_turquoise$estimate, cor_mp17_turquoise$p.value),
  sprintf("MP18 signature\nr = %.2f, p = %.1e", cor_mp18_turquoise$estimate, cor_mp18_turquoise$p.value),
  sprintf("KEGG signature\nr = %.2f, p = %.1e", cor_kegg_turquoise$estimate, cor_kegg_turquoise$p.value),
  sprintf("Assou signature\nr = %.2f, p = %.1e", cor_assou_turquoise$estimate, cor_assou_turquoise$p.value),
  sprintf("Wong signature\nr = %.2f, p = %.1e", cor_wong_turquoise$estimate, cor_wong_turquoise$p.value),
  sprintf("Plurinet signature\nr = %.2f, p = %.1e", cor_plurinet_turquoise$estimate, cor_plurinet_turquoise$p.value),
  sprintf("Benporath signature\nr = %.2f, p = %.1e", cor_benporath_turquoise$estimate, cor_benporath_turquoise$p.value)
)

# ------------------------------------------------------------------
# 6. COLORES Y PLOTS
# ------------------------------------------------------------------

heat_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# 📊 turquoise
pheatmap(mat_turquoise_scaled,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = FALSE,
         show_rownames = TRUE,
         labels_row = labels_turquoise,
         fontsize_row = 9,
         color = heat_colors,
         main = "TURQUOISE MODULE:\n CORRELATION WITH GENETIC SIGNATURES")


## obtenemos los module eigengenes de los modulos de interes 

modulos_interesCOR_m <- c("midnightblue")
MEmodules_AURORA_m <- paste0("ME",modulos_interesCOR_m)

eigen_m <- bwnet$MEs[,MEmodules_AURORA_m]

save(aurora_metastasis, gsva_aurora_met, MEmodules_AURORA_m,expr_midnightblue_m,eigen_m,
     file = "C:/Users/Suruxx/Documents/TFM/data/AURORA_MODULOS_MET.RData")

## -------------------------------------------------------------------
# 7. CORRELACION ENTRE MÓDULOS 
## -------------------------------------------------------------------

##dendograma modulos de interés
plotEigengeneNetworks(eigen_p,
                      names(eigen_p),
                      excludeGrey = TRUE, greyLabel = "grey",
                      plotDendrograms = TRUE, plotHeatmaps = FALSE,
                      colorLabels = TRUE, signed = TRUE,
                      
                      plotAdjacency = FALSE)  # representar correlacion en vez de adyacencia 

## heatmap modulos de interés
plotEigengeneNetworks(eigen_p,
                      names(eigen_p),
                      excludeGrey = TRUE, greyLabel = "grey",
                      plotDendrograms = FALSE, plotHeatmaps = TRUE,
                      colorLabels = TRUE, signed = TRUE,
                      plotAdjacency = FALSE)  # representar correlacion en vez de adyacencia 


## ANÁLISIS ENRIQUECIMIENTO GO 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("enrichplot")
install.packages("ggplot2")

library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

load("C:/Users/Suruxx/Documents/TFM/data/lista_genes_tcga.RData")
load("C:/Users/Suruxx/Documents/TFM/data/lista_genes_aurora_prim.RData")
load("C:/Users/Suruxx/Documents/TFM/data/lista_genes_aurora_met.RData")

## --------------------------------------------------------------------------
# 1. TCGA
## --------------------------------------------------------------------------

## transformamos los simbolos a IDs Entrez
entrez_genes <- bitr(genes_relevantes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
entrez_ids <- entrez_genes$ENTREZID

## análisis de enriquecimiento
ego_BP <- enrichGO(gene          = entrez_ids,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = "BP", 
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2)

ego_CC <- enrichGO(gene          = entrez_ids,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "CC", 
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2)

ego_MF <- enrichGO(gene          = entrez_ids,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "MF", 
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2)


# Barplot
barplot(ego_BP, showCategory = 10) + ggtitle("GO Biological Process Enrichment TCGA data")
barplot(ego_CC, showCategory = 10) + ggtitle("GO Cellular Component Enrichment TCGA data")
barplot(ego_MF, showCategory = 10) + ggtitle("GO Molecular Function Enrichment TCGA data")



## ---------------------------------------------------------------------
# 2. AURORA METASTASIS
## ---------------------------------------------------------------------

### cambiamos algunos genes por su simbolo oficial (consultado en NCBI)
genes_relevantes_met[grep("HIST2H3C",genes_relevantes_met)] <- "H3C14"
genes_relevantes_met[grep("HIST2H3A",genes_relevantes_met)] <- "H3C15"
genes_relevantes_met[grep("HIST1H3B",genes_relevantes_met)] <- "H3C2"
genes_relevantes_met[grep("HIST1H2BL",genes_relevantes_met)] <- "H2BC13"
genes_relevantes_met[grep("HIST1H2AI",genes_relevantes_met)] <- "H2AC13"
genes_relevantes_met[grep("HIST1H3H",genes_relevantes_met)] <- "H3C10"
genes_relevantes_met[grep("HIST1H4J",genes_relevantes_met)] <- "H4C11"
genes_relevantes_met[grep("HIST1H2AL",genes_relevantes_met)] <- "H2AC16"
genes_relevantes_met[grep("HIST1H2BO",genes_relevantes_met)] <- "H2BC17"
genes_relevantes_met[grep("H2AFX",genes_relevantes_met)] <- "H2AX"
genes_relevantes_met[grep("LINC00565",genes_relevantes_met)] <- "SWINGN"

entrez_genes <- bitr(genes_relevantes_met, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
entrez_ids_met <- entrez_genes$ENTREZID

ego_BP_met <- enrichGO(gene          = entrez_ids_met,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "BP", 
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2)

ego_CC_met <- enrichGO(gene          = entrez_ids_met,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "CC", 
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2)

ego_MF_met <- enrichGO(gene          = entrez_ids_met,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENTREZID',
                   ont           = "MF", 
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.2)

# Barplot
barplot(ego_MF_met, showCategory = 15) + ggtitle("GO Molecular Function Enrichment")



## ----------------------------------------------------------------------------
# 3. AURORA PRIMARIOS
## ----------------------------------------------------------------------------

entrez_genes <- bitr(genes_relevantes_prim, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
entrez_ids_prim <- entrez_genes$ENTREZID

ego_BP_prim <- enrichGO(gene          = entrez_ids_prim,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "BP", 
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.2)

ego_CC_prim <- enrichGO(gene          = entrez_ids_prim,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "CC", 
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.2)

ego_MF_prim <- enrichGO(gene          = entrez_ids_prim,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "MF", 
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.2)

# Barplot
barplot(ego_BP_met, showCategory = 10) + ggtitle("GO Biological Process Enrichment AURORA metastasis")
barplot(ego_CC_met, showCategory = 10) + ggtitle("GO Cellular Component Enrichment AURORA metastasis")
barplot(ego_MF_met, showCategory = 10) + ggtitle("GO Molecular Function Enrichment AURORA metastasis")



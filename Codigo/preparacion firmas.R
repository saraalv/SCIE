## FIRMAS GENICAS

# Cargamos los datos 
firmas <- read.csv2("C:/Users/Suruxx/Documents/TFM/data/firmasgenicas.csv")

# Separamos cada una de las firmas en un objeto 

ISDS_signature <- firmas$ISDS_signature[1:27]
MP17 <- firmas$MP17.Interferon.MHC.II..I.[1:50]
MP18 <- firmas$MP18.Interferon.MHC.II..II.[1:50]
KEGG <- firmas$KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION[1:88]
assou <- firmas$Hs_ESC_Assou
wong <- firmas$Hs_ESC_Wong[1:335]
plurinet <- firmas$SC_Plurinet_computational[1:299]
benporath <- firmas$ES_exp1_benporath[1:380]

# creamos una lista que contenga las 8 firmas génicas 

firmas <- list(ISDS_signature,MP17,MP18,KEGG,assou,wong,plurinet,benporath)
names(firmas) <- c("ISDS","MP17","MP18","KEGG","assou","wong","plurinet","benporath")

#cargamos los metadatos

firmas_met <- read.csv2("C:/Users/Suruxx/Documents/TFM/data/firmasgenicas_met.csv")

#guardamos la lista de firmas y sus metadatos en un archivo RData

save(firmas, firmas_met, file = "C:/Users/Suruxx/Documents/TFM/data/firmasgenicas.RData")

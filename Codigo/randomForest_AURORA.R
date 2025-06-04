
## Feature Selection en dataset AURORA

library(randomForest)
library(caret)

load(file = "C:/Users/Suruxx/Documents/TFM/data/aurora_dataset.RData")
load("C:/Users/Suruxx/Documents/TFM/data/AURORA_MODULOS_MET.RData")
load("C:/Users/Suruxx/Documents/TFM/data/AURORA_MODULOS.RData")
load("C:/Users/Suruxx/Documents/TFM/data/AURORA_MODULOS_PRIM.RData")
load("C:/Users/Suruxx/Documents/TFM/data/df_modulos_aurora_prim.RData")
load("C:/Users/Suruxx/Documents/TFM/data/df_modulos_aurora.RData")


##---------------------------------------------------------------------
# 1. División de las muestras según scores
##---------------------------------------------------------------------

scores_mat <- cbind(
  stemness = gsva_aurora_met["benporath",],
  immunogenicity = gsva_aurora_met["MP17",])

rownames(scores_mat) <- colnames(aurora_metastasis)
plot(scores_mat)

high_stem <- c()
low_stem <- c()
high_immune <- c()
low_immune <- c()


## inicialmente se iba a identificar como high stem/immune las muestras con gsva 
## mayor de el tercer cuantil y como low stem/immune las que tengan el gsva por 
## debajo del tercer cuantil, como en el caso de AURORA nos quedariamos con una n
## muy baja (11 para primarios y 20 para metástasis), marcaremos el umbral en la mediana.

for (i in colnames(aurora_metasasis)) {
  if (gsva_aurora_met["benporath",i] >= median(gsva_aurora_met["benporath",])){
    high_stem <- c(high_stem, i)
    
  }
  if (gsva_aurora_met["benporath",i] <= median(gsva_aurora_met["benporath",])) {
    low_stem <- c(low_stem, i)
    
    
  }
  if (gsva_aurora_met["MP17",i] >= median(gsva_aurora_met["MP17",])){
    high_immune <- c(high_immune,i)
    
  }
  if (gsva_aurora_met["MP17", i] <= median(gsva_aurora_met["MP17",])){
    low_immune <- c(low_immune,i)
    
  }
  
}

## asignamos a cada muestra una clase 

met_ids <- c()
clase <- c()
for (i in colnames(aurora_metastasis)) {
  if (i %in% high_stem[high_stem %in% high_immune]) {
    clase <- c(clase,"HSHI")
    met_ids <- c(met_ids,i)
  }
  if (i %in% high_stem[high_stem %in% low_immune]){
    clase <- c(clase, "HSLI")
    met_ids <- c(met_ids,i)
  }
  if (i %in% low_stem[low_stem %in% high_immune]){
    clase <- c(clase,"LSHI")
    met_ids <- c(met_ids,i)
  }
  if (i %in% low_stem[low_stem %in% low_immune]){
    clase <- c(clase,"LSLI")
    met_ids <- c(met_ids,i)
  }
}

## creamos un dataframe con los datos de expresion de todos los genes de interés
## y el tipo de tumor al que corresponden

###Primarios
modulos_interes <- gsub("ME","",MEmodules_AURORA_p)
genes_id <- subset(modules_df_prim$gene_id, modules_df_prim$module %in% modulos_interes)

###Metastasis
genes_id <- rownames(expr_midnightblue_m)

clas_df <- t(aurora_data[genes_id,met_ids])
clas_df <- as.data.frame(clas_df)
clas_df <- cbind(clas_df,clase)
clas_df$clase <- factor(clas_df$clase)


### Renombramos algunas variables que dan error en Random Forest: 

colnames(clas_df) <- gsub("-","_",colnames(clas_df))

##----------------------------------------------------------------------
# 2 . Creación del modelos de Random Forest
##---------------------------------------------------------------------

## dividimos en training y test data set 

set.seed(1412)
ntrain <- (2/3)*length(clas_df$clase)
numtrain <- sample(length(clas_df$clase),ntrain)
au_train <- clas_df[numtrain,]
au_test <- clas_df[-numtrain,]
labels_train <- clas_df[numtrain,"clase"]
labels_test <- clas_df[-numtrain,"clase"]

## creamos el modelo de Random Forest

control <- trainControl(method = "cv", number = 10, verboseIter = TRUE)

rf_met <- train( clase ~.,
                  data = au_train,
                  method = "rf",
                  trControl = control,
                  tuneLength = 5)


saveRDS(rf_met, file ="C:/Users/Suruxx/Documents/TFM/data/random_forest_met.rds" )

## -------------------------------------------------------------------------
# 3. Evaluamos su actuación como clasificador  
## -------------------------------------------------------------------------

# evaluamos su rendimiento clasificando nuevos datos 

rf_pred <- predict(rf_met, au_test)
cm <- confusionMatrix(rf_pred,labels_test)
cm_df <- as.data.frame(cm$table)

# representamo la matriz de confusion 

ggplot(data = cm_df, aes(x = Prediction, y = Reference, fill = Freq)) +
  geom_tile(color = "honeydew") +
  geom_text(aes(label = Freq), size = 4) +
  scale_fill_gradient(low = "honeydew", high = "steelblue") +
  labs(title = "Matriz de Confusión AURORA Metástasis",
       x = "Predicción",
       y = "Valor Real") +
  theme_minimal() +
  theme(legend.position = "none")

## ----------------------------------------------------------------
# 4. SELECCIÓN DE GENES 
## ----------------------------------------------------------------

# seleccionamos los genes con una importancia mayor de 45 

imp <- varImp(rf_met)$importance
genes_relevantes_met <- rownames(imp)[imp$Overall >= 45]
length(genes_relevantes_met) 

# entrenamos un nuevo modelo con los genes seleccionados 

clas_sel <- clas_df[,colnames(clas_df) %in% c(genes_relevantes_met,"clase")]
au_train_sel <- clas_sel[numtrain,]
au_test_sel <- clas_sel[-numtrain,]
labels_train_s <- clas_sel[numtrain,"clase"]
labels_test_s <- clas_sel[-numtrain,"clase"]

rf_rel <- train( clase ~.,
                 data = au_train_sel,
                 method = "rf",
                 trControl = control,
                 tuneLength = 5)

# evaluamos su rendimiento

pred <- predict(rf_rel, au_test_sel)
cm_s <- confusionMatrix(pred, labels_test_s)
cm_df_s <- as.data.frame(cm_s$table)

ggplot(data = cm_df_s, aes(x = Prediction, y = Reference, fill = Freq)) +
  geom_tile(color = "honeydew") +
  geom_text(aes(label = Freq), size = 4) +
  scale_fill_gradient(low = "honeydew", high = "steelblue") +
  labs(title = "Matriz de Confusión AURORA Genes Seleccionados",
       x = "Predicción",
       y = "Valor Real") +
  theme_minimal() +
  theme(legend.position = "none")

genes_relevantes_au <- gsub("_","-",genes_relevantes_met) # volvemos los IDs a su formato original 

# guardamos la lista de genes mas informativos 

save(genes_relevantes_met, file = "C:/Users/Suruxx/Documents/TFM/data/lista_genes_aurora_met.RData")

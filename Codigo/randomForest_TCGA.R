
## SELECCIÓN DE GENES MÁS INFORMATIVOS 

library(randomForest)
library(caret)

load(file = "C:/Users/Suruxx/Documents/TFM/data/tcga_dataset.RData")
load(file = "C:/Users/Suruxx/Documents/TFM/data/TCGA_MODULOS.RData")
load("C:/Users/Suruxx/Documents/TFM/data/modulos_df_tcga.RData")

load("C:/Users/Suruxx/Documents/TFM/data/firmasgenicas.RData")

##---------------------------------------------------------------------
# 1. División de las muestras según scores
##---------------------------------------------------------------------

scores_mat <- cbind(
  stemness = gsva_tcga["benporath",],
  immunogenicity = gsva_tcga["MP17",])

rownames(scores_mat) <- colnames(tcga_data)
plot(scores_mat)

high_stem <- c()
low_stem <- c()
high_immune <- c()
low_immune <- c()


## identificamos como high stem/immune las muestras con gsva 
## mayor de el tercer cuantil y como low stem/immune las que
## tengan el gsva por debajo del tercer cuantil.

for (i in colnames(tcga_data)) {
  if (gsva_tcga["benporath",i] >= quantile(gsva_tcga["benporath",], probs = 0.75)){
    high_stem <- c(high_stem, i)
    
  }
  if (gsva_tcga["benporath",i] <= quantile(gsva_tcga["benporath",], probs = 0.25)) {
    low_stem <- c(low_stem, i)
    
    
  }
  if (gsva_tcga["MP17",i] >= quantile(gsva_tcga["MP17",], probs = 0.75)){
    high_immune <- c(high_immune,i)
    
  }
  if (gsva_tcga["MP17", i] <= quantile(gsva_tcga["MP17",], probs = 0.25)){
    low_immune <- c(low_immune,i)
    
  }
  
}

## asignamos a cada muestra una clase 

tcga_ids <- c()
clase <- c()
for (i in colnames(tcga_data)) {
  if (i %in% high_stem[high_stem %in% high_immune]) {
    clase <- c(clase,"HSHI")
    tcga_ids <- c(tcga_ids,i)
  }
  if (i %in% high_stem[high_stem %in% low_immune]){
    clase <- c(clase, "HSLI")
    tcga_ids <- c(tcga_ids,i)
  }
  if (i %in% low_stem[low_stem %in% high_immune]){
    clase <- c(clase,"LSHI")
    tcga_ids <- c(tcga_ids,i)
  }
  if (i %in% low_stem[low_stem %in% low_immune]){
    clase <- c(clase,"LSLI")
    tcga_ids <- c(tcga_ids,i)
  }
}

## creamos un dataframe con los datos de expresion de todos los genes de interés
## y el tipo de tumor al que corresponden

modulos_interes <- gsub("ME","",MEmodules_tcga)
genes_id <- subset(modules_df_tcga$gene_id, modules_df_tcga$module %in% modulos_interes)

clas_df <- t(tcga_data[genes_id,tcga_ids])
clas_df <- as.data.frame(clas_df)
clas_df <- cbind(clas_df,clase)
clas_df$clase <- factor(clas_df$clase)


### Renombramos algunas variables que dan error en Random Forest: 

colnames(clas_df) <- gsub("-","_",colnames(clas_df))

##----------------------------------------------------------------------
# 2. Creación del modelo de Random Forest
##---------------------------------------------------------------------

## dividimos en training y test data set 

set.seed(1412)
ntrain <- (2/3)*length(clas_df$clase)
numtrain <- sample(length(clas_df$clase),ntrain)
tcga_train <- clas_df[numtrain,]
tcga_test <- clas_df[-numtrain,]
labels_train <- clas_df[numtrain,"clase"]
labels_test <- clas_df[-numtrain,"clase"]

## creamos el modelo de Random Forest

control <- trainControl(method = "cv", number = 10, verboseIter = TRUE)

rf_tcga <- train( clase ~.,
             data = tcga_train,
             method = "rf",
             trControl = control,
             tuneLength = 5)


saveRDS(rf_tcga, file ="C:/Users/Suruxx/Documents/TFM/data/random_forest_tcga.rds" )

## -------------------------------------------------------------------
# 3. Evaluamos su actuación como clasificador  
## -------------------------------------------------------------------

# evaluamos su rendimiento clasificando nuevos datos 

rf_pred <- predict(rf_tcga, tcga_test)
cm <- confusionMatrix(rf_pred,labels_test)
cm_df <- as.data.frame(cm$table)

# representamos graficamente la matriz de confusion 

ggplot(data = cm_df, aes(x = Prediction, y = Reference, fill = Freq)) +
  geom_tile(color = "honeydew") +
  geom_text(aes(label = Freq), size = 4) +
  scale_fill_gradient(low = "honeydew", high = "steelblue") +
  labs(title = "Matriz de Confusión TCGA",
       x = "Predicción",
       y = "Valor Real") +
  theme_minimal() +
  theme(legend.position = "none")

## ----------------------------------------------------------------
# 4. SELECCIÓN DE GENES 
## ----------------------------------------------------------------

# seleccionamos los genes con una importancia mayor de 45 

imp <- varImp(rf_tcga)$importance
genes_relevantes <- rownames(imp)[imp$Overall >= 45]
length(genes_relevantes) 

# entrenamos un nuevo modelo con los genes seleccionados 

clas_sel <- clas_df[,colnames(clas_df) %in% c(genes_relevantes,"clase")]

tcga_train_sel <- clas_sel[numtrain,]
tcga_test_sel <- clas_sel[-numtrain,]
labels_train_s <- clas_sel[numtrain,"clase"]
labels_test_s <- clas_sel[-numtrain,"clase"]

rf_rel <- randomForest(clase ~., tcga_train_sel)

# evaluamos su rendimiento
pred <- predict(rf_rel, tcga_test_sel)
cm_s <- confusionMatrix(pred, labels_test_s)
cm_df_s <- as.data.frame(cm_s$table)

ggplot(data = cm_df_s, aes(x = Prediction, y = Reference, fill = Freq)) +
  geom_tile(color = "honeydew") +
  geom_text(aes(label = Freq), size = 4) +
  scale_fill_gradient(low = "honeydew", high = "steelblue") +
  labs(title = "Matriz de Confusión TCGA Genes Seleccionados",
       x = "Predicción",
       y = "Valor Real") +
  theme_minimal() +
  theme(legend.position = "none")

genes_relevantes <- gsub("_","-",genes_relevantes) # volvemos los IDs a su formato original 

# guardamos la lista de genes mas informativos 

save(genes_relevantes, file = "C:/Users/Suruxx/Documents/TFM/data/lista_genes_tcga.RData")


### Author Jacob Schladebach, 30/10/2024

library(lidR)
library(dplyr)
require(caret)
require(MLmetrics)


training <- read.csv2("training_data.csv")

equal_classes <- training %>%
  group_by(Species) %>%
  mutate(row_id = row_number()) %>%  
  filter(row_id <= 50) %>%            
  sample_n(size = min(n(), 50)) %>%   
  ungroup() %>%
  select(-row_id)  

train_set<-data.frame(Zq999 = equal_classes$Zq999, NGRDI = equal_classes$NGRDI, VARI = equal_classes$VARI, Zq99= equal_classes$Zq99, vol_concave = equal_classes$vol_concave,
                      density_concave = equal_classes$density_concave, CV_Z = equal_classes$CV_Z, CRR=equal_classes$CRR,  vertical_variability=equal_classes$vertical_variability,
                      pointedness= equal_classes$pointedness, top_angle = equal_classes$top_angle, widest_distance=equal_classes$widest_distance,
                      relative_height_widest=equal_classes$relative_height_widest, Species= equal_classes$Species)

train_set$Species <- as.factor(train_set$Species)
train_set<-na.omit(train_set)

train_set$Species<-gsub("3","PIMA",train_set$Species)
train_set$Species<-gsub("4","PIGL",train_set$Species)
train_set$Species<-gsub("5","PICO",train_set$Species)
train_set$Species<-gsub("6","ABLA",train_set$Species)
train_set$Species<-gsub("7","LALA",train_set$Species)
train_set$Species<-gsub("8","BEPA",train_set$Species)
train_set$Species<-gsub("9","POTR",train_set$Species)
train_set$Species<-gsub("10","POBA",train_set$Species)
train_set$Species<-gsub("11","PISI",train_set$Species)
train_set$Species<-gsub("12","ALGL",train_set$Species)
train_set$Species<-gsub("13","BENE",train_set$Species)

train_control <- trainControl(method = "cv", number = 10, 
                              classProbs = TRUE, 
                              summaryFunction = multiClassSummary)

rf_model_kfold<- train(Species ~ Zq999 + Zq99  + vol_concave + VARI + NGRDI +density_concave + CV_Z + CRR + pointedness+ vertical_variability + top_angle+ widest_distance  + relative_height_widest,
                        data = train_set,
                        ntree = 500,
                        method = "rf", trControl = train_control)
print(rf_model_kfold)

saveRDS(rf_model_kfold, "Colour_model_NA.rds")

importance <- varImp(rf_model_kfold, scale = FALSE)
plot(importance)

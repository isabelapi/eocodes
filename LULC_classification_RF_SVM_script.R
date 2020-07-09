###################
#Script to make Random Forest (RF) and Suport Vector Machine(SVM) 
#on Land-use/Land-cover (LULC) classification on remote sensing images

###LULC calssification of HYSPEX###
{
  require(raster)
  require(caret)
  require(rgdal)
  require(randomForest)
  require(caTools)
  require(e1071)
  setwd("F:/MSc._Environmental_Science/MSc_ES_Semester_2/Adv_Remote_Sensing/Term_paper_Adv_Remote_Sensing/Part_2-Land Use Classification/Classification")

  ## RF ##
  pca.image <- brick("E:/MSc._Environmental_Science/MSc_ES_Semester_2/Adv_Remote_Sensing/Term_paper_Adv_Remote_Sensing/Part_2-Land Use Classification/Classification/PCA_Cov_19b")
  plotRGB(pca.image, 3, 2 ,1, stretch = "lin")
  training <- readOGR("E:/MSc._Environmental_Science/MSc_ES_Semester_2/Adv_Remote_Sensing/Term_paper_Adv_Remote_Sensing/Part_2-Land Use Classification/Classification", "roi_merge")
  
  #Convert levels to factor, convert classes from text to numbers
  training$Class <- as.factor(training$Class)
  Class = levels(training$Class)
  for (i in 1:length(Class)){
    levels(training$Class)[levels(training$Class)==Class[[i]]] <- i
  }
  
  #extract training 
  reflectance <- extract(pca.image, training, df = T)
  reflectance$Class <- as.factor(training$Class[reflectance$ID])
  colnames(reflectance) <- c("ID", "pc.1", "pc.2", "pc.3", "pc.4", "pc.5", "pc.6", "pc.7", "pc.8", "pc.9"
                             , "pc.10", "pc.11", "pc.12", "pc.13", "pc.14", "pc.15", "pc.16", "pc.17", "pc.18", "pc.19","Class")
  
  #split data into training (40%), validation (60%)
  set.seed(1000)
  sample <- sample.split(reflectance, SplitRatio = 0.3)
  train <- subset (reflectance, sample == TRUE)
  validation <- subset (reflectance, sample == FALSE) 
  validation <- validation[,-1] # remove ID column
  nrow(train)/nrow(validation)
  train <- train[, -1] # remove ID column
  
  #Grid searches RF parameters trees and variable split
  ntrees<-seq(100,1000,by=100)
  mtrys<-seq(2,6,by=1)
  kappamatrix<-matrix(0,nrow=length(ntrees),ncol=length(mtrys))
  rownames(kappamatrix)<-ntrees
  colnames(kappamatrix)<-mtrys
  for (i in 1:length(ntrees)){
    for (j in 1:length(mtrys)){
      mtry<-mtrys[j]
      ntree<-ntrees[i]
      cat(sprintf("ntree: %f mtry:%f",ntree,mtry))
      model<-randomForest(Class ~., data = train, importance = T, ntree = ntree, mtry = mtry)
      kappa<-confusionMatrix(validation$Class,predict(model,validation))$overall[2]
      kappamatrix[as.character(ntree),as.character(mtry)]<-kappa
      cat(sprintf(" kappa: %f \n",kappa))
    }
  }
  plot(kappamatrix)
  
  #get maximum kappa
  k <- arrayInd(which.max(kappamatrix), dim(kappamatrix))
  ntree = as.numeric(rownames(kappamatrix)[k[,1]])
  mtry = as.numeric(colnames(kappamatrix)[k[,2]])
  
  beginCluster()
  
  #train model
  rf <- randomForest(Class ~., data = train, importance = T, mtry = mtry, ntree = ntree)
  print(rf)
  class <- predict(rf, newdata = validation)
  
  #accuracy assessment
  cf <- confusionMatrix(class, validation$Class) 
  print(cf)
  AccAss_rf <- t(cf$table) #use for AccAss
  
  #convert to data frame
  for (i in 1:dim(pca.image)[3]){
    assign(paste("pc.", i, sep = ""), as.vector(as.matrix(pca.image[[i]])))
  }
  data <- data.frame(pc.1, pc.2, pc.3, pc.4, pc.5, pc.6, pc.7, pc.8, pc.9
                     , pc.10, pc.11, pc.12, pc.13, pc.14, pc.15, pc.16, pc.17, pc.18, pc.19)
  
  #final RF classification 
  final.class <- predict(rf,newdata=data)
  pred_num<-as.numeric(levels(final.class)[final.class])
  table_pix_rf <- table(pred_num) # number of pixels per class
  
  #convert final classfication to matrix - raster
  pred_matrix<-matrix(pred_num,ncol=2000,nrow=2000)
  pred_raster<-raster(pred_matrix)
  
  #plot RF classification raster
  plot(pred_raster)
  crs(pred_raster)<-crs(pca.image)
  extent(pred_raster)<-extent(pca.image)

  ## SVM ##
  
  #Grid search for optimal Radial SVM parameters cost and gamma
  costs<-seq(0.2,2,by=0.2)
  gammas<-seq(0.01,0.1,by=0.01)
  kappamatrix<-matrix(0,nrow=length(costs),ncol=length(gammas))
  rownames(kappamatrix)<-costs
  colnames(kappamatrix)<-gammas
  for (i in 1:length(costs)){
    for (j in 1:length(gammas)){
      gamma<-gammas[j]
      cost<-costs[i]
      cat(sprintf("cost: %f gamma:%f",cost,gamma))
      model<-svm(Class~.,data=train,kernel="radial",gamma=gamma,cost=cost)
      kappa<-confusionMatrix(validation$Class,predict(model,validation))$overall[2]
      kappamatrix[as.character(cost),as.character(gamma)]<-kappa
      cat(sprintf(" kappa: %f \n",kappa))
    }
  }
  plot(kappamatrix)
  
  #get maximum kappa
  k <- arrayInd(which.max(kappamatrix), dim(kappamatrix))
  cost = as.numeric(rownames(kappamatrix)[k[,1]])
  gamma = as.numeric(colnames(kappamatrix)[k[,2]])
  
  #train model
  svm <- svm (Class ~ . , data = train, kernel = "radial", gamma = gamma, cost = cost)
  
  #accuracy assessment
  predict_svm <- predict (svm, validation)
  cf_svm <- confusionMatrix (predict_svm, validation$Class)
  print(cf_svm)
  AccAss_svm <- t(cf_svm$table) #use for AccAss
  
  
  #SVM land use classification
  ClassSVM <- predict(model,newdata=data)
  SVM_num<-as.numeric(levels(ClassSVM)[ClassSVM])
  table_pix_svm <- table(SVM_num) # number of pixels per class
  
  #convert final classfication to matrix - raster
  SVM_matrix<-matrix(SVM_num,ncol=2000,nrow=2000)
  SVM_raster<-raster(SVM_matrix)
  
  #plot SVM classification raster
  plot(SVM_raster)
  crs(SVM_raster)<-crs(pca.image)
  extent(SVM_raster)<-extent(pca.image)
  
  #export raster and Accuracy Assesment of both classifications of HySpex image
  writeRaster(pred_raster, filename = "rf_class_hyspex_1", format = "GTiff", overwrite = TRUE)
  write.table(AccAss_rf, file = "rf_hyspex_AccAss.csv", sep = " ", row.names = TRUE, col.names = TRUE)

  writeRaster(SVM_raster, filename = "SVM_class_hyspex_1", format = "GTiff", overwrite = TRUE)
  write.table(AccAss_svm, file = "svm_hyspex_AccAss.csv", sep = " ", row.names = TRUE, col.names = TRUE)
  
  endCluster()
}


###LULC calssification of Landsat 8 multitemporal stack###

{
  require(raster)
  require(caret)
  require(rgdal)
  require(randomForest)
  require(caTools)
  require(e1071)
  
  land.image <- brick("F:/MSc._Environmental_Science/MSc_ES_Semester_2/Adv_Remote_Sensing/Term_paper_Adv_Remote_Sensing/Part_2-Land Use Classification/Classification/Landsat classification/landsat_stack_PCA_spectralsubset.tif")
  plotRGB(land.image, 3, 2 ,1, stretch = "lin")
  training <- readOGR("F:/MSc._Environmental_Science/MSc_ES_Semester_2/Adv_Remote_Sensing/Term_paper_Adv_Remote_Sensing/Part_2-Land Use Classification/Classification", "roi_merge")
  
  #Convert levels to factor, convert classes from text to numbers
  training$Class <- as.factor(training$Class)
  Class = levels(training$Class)
  
  for (i in 1:length(Class)){
    levels(training$Class)[levels(training$Class)==Class[[i]]] <- i
  }
  
  #extract training 
  reflectance <- extract(land.image, training, df = T)
  reflectance$Class <- as.factor(training$Class[reflectance$ID])
  colnames(reflectance) <- c("ID", "b1", "b2", "b3", "b4", "b5", "Class")
  
  #split data into training (40%), validation (60%)
  set.seed(1000)
  reflectance$ID <- 1:nrow (reflectance)
  Classes <- as.vector(as.factor(2:11))
  temp <- reflectance[reflectance$Class == "1",]
  tempID <- sample.split(temp$Class, SplitRatio = 0.3)
  tempSample <- temp[tempID,]
  train <- tempSample
  
  for (i in Classes) {
    temp <- reflectance[reflectance$Class == i,]
    tempID <- sample.split(temp$Class, SplitRatio = 0.3)
    tempSample <- temp[tempID,]
    train <- rbind (train, tempSample)
  }
  
  validationID <- reflectance$ID [-train$ID]
  validation <- reflectance [validationID, -1] #remove ID column
  nrow(train)/nrow(validation) #to verify splitting 
  train <- train[, -1] # remove ID column
  
  #Grid searches RF parameters trees and variable split
  ntrees<-seq(100,1000,by=100)
  mtrys<-seq(2,6,by=1)
  kappamatrix<-matrix(0,nrow=length(ntrees),ncol=length(mtrys))
  rownames(kappamatrix)<-ntrees
  colnames(kappamatrix)<-mtrys
  for (i in 1:length(ntrees)){
    for (j in 1:length(mtrys)){
      mtry<-mtrys[j]
      ntree<-ntrees[i]
      cat(sprintf("ntree: %f mtry:%f",ntree,mtry))
      model<-randomForest(Class ~., data = train, importance = T, ntree = ntree, mtry = mtry)
      kappa<-confusionMatrix(validation$Class,predict(model,validation))$overall[2]
      kappamatrix[as.character(ntree),as.character(mtry)]<-kappa
      cat(sprintf(" kappa: %f \n",kappa))
    }
  }
  plot(kappamatrix)
  
  #get maximum kappa
  k <- arrayInd(which.max(kappamatrix), dim(kappamatrix))
  ntree = as.numeric(rownames(kappamatrix)[k[,1]])
  mtry = as.numeric(colnames(kappamatrix)[k[,2]])
  
  beginCluster()
  
  #train model
  rf <- randomForest(Class ~., data = train, importance = T, mtry = mtry, ntree = ntree)
  print(rf)
  class <- predict(rf, newdata = validation)
  
  #accuracy assessment
  cf <- confusionMatrix(class, validation$Class) #######
  print(cf)
  AccAss <- t(cf$table) #use for AccAss
  
  #convert to data frame
  for (i in 1:dim(land.image)[3]){
    assign(paste("land.", i, sep = ""), as.vector(as.matrix(land.image[[i]])))
  }
  
  data <- data.frame(land.1, land.2, land.3, land.4, land.5)
  colnames(data) <- colnames (reflectance[,2:6])
  
  #final RF classification 
  final.class <- predict(rf,newdata=data)
  pred_num<-as.numeric(levels(final.class)[final.class])
  table_pix_rf <- table(pred_num) # number of pixels per class
  
  #convert final classfication to matrix - raster
  pred_matrix<-matrix(pred_num,ncol=1000,nrow=1000) #this should be same as dimention of ncols and nrows of original image
  pred_raster<-raster(pred_matrix)
  
  #plot classification raster
  plot(pred_raster)
  crs(pred_raster)<-crs(land.image)
  extent(pred_raster)<-extent(land.image)
  
  #export raster and Accuracy Assessment of RF classification
  writeRaster(pred_raster, filename = "rf_class_landsat", format = "GTiff", overwrite = TRUE)
  write.table(AccAss, file = "rf_class_landsat_AccAss.csv", sep = " ", row.names = TRUE, col.names = TRUE)
  
  endCluster()
}
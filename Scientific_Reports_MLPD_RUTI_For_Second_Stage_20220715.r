# The code for the paper
# Machine learning to predict the development of recurrent urinary tract infection related to single uropathogen, Escherichia coli
#
#
#Shuen-Lin Jeng1, Zi-Jing Huang2, Deng-Chi Yang3, Ching-Hao Teng4*, Ming-Cheng Wang5*
#
#1Department of Statistics, Institute of Data Science, and Center for Innovative FinTech Business Models, National Cheng Kung University, Tainan, Taiwan
#2 Department of Statistics, National Cheng Kung University, Tainan, Taiwan
#3Department of Geriatrics and Gerontology, National Cheng Kung University Hospital, College of Medicine, National Cheng Kung University, Tainan, Taiwan.
#4Institute of Molecular Medicine, College of Medicine, National Cheng Kung University, Tainan, Taiwan; Institute of Basic Medical Sciences, College of Medicine, National Cheng Kung University, Tainan, Taiwan; Center of Infectious Disease and Signaling Research, National Cheng Kung University, Tainan, Taiwan.
#5Division of Nephrology, Department of Internal Medicine, National Cheng Kung University Hospital, College of Medicine, National Cheng Kung University, Tainan, Taiwan.
#

library(readxl)
library(magrittr)
library(randomForest)
library(rpart)
library(rpart.plot)

#===========================================================
# Second stage Analysis

############################################################
#####             Read data file and data processing   #####
############################################################
UTIdata.df <- read_xlsx("C:/Users/88692/Desktop/Scientific_Reports_MLPD_RUTI/Original_RUTI_Data_Without_PatientInfo.xlsx",sheet=1)
UTIdata.df %<>% as.data.frame(.)
colname.v <- c("Bacterial_Name",
               paste0("Anti",1:37),
               "Ecoli_No",
               paste0("Gene",1:33), 
               "Recurrent",
               "Gender","Age",
               "Hospital_day","Dead",
               "hospital_long","Clinical_features",
               "UTI_pos",
               paste0("Dis",1:12),
               "Pre_hos_2y",
               "Pre_UTI_ER_2y","Pre_UTI_hos_2y",
               "Symptomatic",
               paste0("UTI_symptom",1:8), 
               "BloodWBC","Creatinine","UWBC","UBact","URBC","ER_H_OPD") 
colnames(UTIdata.df) <- colname.v

UTIdata.df$Age %<>% ifelse(.==767,NA,.)
UTIdata.df$Dead %<>% ifelse(.==9, NA, .) %>% as.factor()
UTIdata.df$Gender %<>% as.factor(.)
UTIdata.df$UTI_pos %<>% ifelse(. %in% c("X","x","0"), NA, .) %>% as.factor()
UTIdata.df$Pre_hos_2y %<>% ifelse(. %in% c("x"), NA, .) %>% as.numeric()
UTIdata.df$Pre_UTI_ER_2y %<>% ifelse(. %in% c("x"), NA, .) %>% as.numeric()
UTIdata.df$Pre_UTI_hos_2y %<>% ifelse(. %in% c("x"), NA, .) %>% as.numeric()
UTIdata.df$UTI_symptom1 %<>% ifelse(. %in% c(2,11), NA, .) 
UTIdata.df$UTI_symptom5 %<>% ifelse(. %in% c("x"), NA, .)
UTIdata.df$UTI_symptom8 %<>% ifelse(. == 19.9, NA, .)
for (i in c(97:104)) {
  UTIdata.df[,i] <- factor(UTIdata.df[,i], levels = c("0","1"))
}
UTIdata.df$BloodWBC %<>% ifelse(. %in% c("x"), NA, .) %>% as.numeric()
UTIdata.df$Creatinine %<>% ifelse(. %in% c("x"), NA, .) %>% as.numeric()
UTIdata.df$UBact %<>% ifelse(. %in% c("x"), NA, .) %>% as.factor()
UTIdata.df$UBact %<>% ordered(., levels = c(0,1,2,3,4))

RBC_level.f <- function(x,dat){
  i=dat[x]
  if(!is.na(i)){
    if(i==0){
      idx=0
    }else if(i>=1 && i<=10){
      idx=1
    }else if(i>=11 && i<=100){
      idx=2
    }else if(i>=101 && i<=1000){
      idx=3
    }else if(i>=1001){
      idx=4
    }else{
      idx=NA
    }
  }else{
    idx=NA
  }
  return(idx)
}

WBC_level.f <- function(x,dat){
  i=dat[x]
  if(!is.na(i)){
    if(i>=0 && i<=10){
      idx=1
    }else if(i>=11 && i<=100){
      idx=2
    }else if(i>=101 && i<=1000){
      idx=3
    }else if(i>=1001){
      idx=4
    }else{
      idx=NA
    }
  }else{
    idx=NA
  }
  return(idx)
}

UTIdata.df$UWBC %<>% ifelse(. %in% c("x","3+"), NA, .)
UTIdata.df$UWBC %<>% ifelse(. %in% c("numerous"), 1001, .)
UTIdata.df$UWBC %<>% ifelse(. %in% c("100+"), 101, .)
UTIdata.df$UWBC %<>% ifelse(. %in% c("1000+"), 1001, .)
UTIdata.df$UWBC %<>% as.numeric(.)
UTIdata.df$UWBC_level <- sapply(1:length(UTIdata.df$UWBC), WBC_level.f, dat=UTIdata.df$UWBC) %>% factor()
UTIdata.df$UWBC_level %<>%  ordered(., levels = c(1,2,3,4))

UTIdata.df$URBC %<>% ifelse(. %in% c("x"), NA, .)
UTIdata.df$URBC %<>% ifelse(. %in% c("numerous"), 1001, .)
UTIdata.df$URBC %<>% ifelse(. %in% c("100+"), 101, .)
UTIdata.df$URBC %<>% ifelse(. %in% c("1000+"), 1001, .)
UTIdata.df$URBC %<>% as.numeric(.)
UTIdata.df$URBC_level <- sapply(1:length(UTIdata.df$URBC), RBC_level.f, dat=UTIdata.df$URBC) %>% factor()
UTIdata.df$URBC_level %<>% ordered(., levels = c(0,1,2,3,4))

UTIdata.df$Bacterial_Name %<>% ifelse(. %in% c("Klebsiella pneumoniae"), NA, .) 
UTIdata.df$Bacterial_Name %<>% ifelse(. %in% c("Escherichia coli"), 1, .) 
UTIdata.df$Bacterial_Name %<>% ifelse(. %in% c("Escherichia coli with ESBL"), 2, .)
UTIdata.df$Bacterial_Name %<>% factor(.)
UTIdata.df$Recurrent %<>% as.factor(.)

## Anti1~37
anti.f <- function(dat){
  dat <- dat %>% as.matrix %>% as.character
  dat <- ifelse(dat=="S", 1, dat)
  dat <- ifelse(dat=="I", 2, dat)
  dat <- ifelse(dat=="R", 3, dat)
  dat <- ifelse(is.na(dat), NA, dat)
  dat <- dat %>% factor(., levels = c("1","2","3"))
  dat %<>% ordered(., levels = c(1,2,3))
  dat %>% as.matrix() %>% return()
}
colnames(UTIdata.df)[2:38]
UTIdata.df[,2:38] <- apply(UTIdata.df[,2:38], 2,anti.f) 
for (i in 2:38) {
  UTIdata.df[,i] <- factor(UTIdata.df[,i], levels = c("1","2","3"))
}

## Gene1~16,18~33
colnames(UTIdata.df)[c(40:55,57:72)]
for (i in c(40:55,57:72)) {
  UTIdata.df[,i] <- factor(UTIdata.df[,i], levels=c("0","1"))
}
## Gene17
colnames(UTIdata.df)[56]
gene17.f <- function(dat){
  dat <- dat %>% as.matrix %>% as.character
  dat <- ifelse(dat=="A", 1, dat)
  dat <- ifelse(dat=="B1", 2, dat)
  dat <- ifelse(dat=="B2", 3, dat)
  dat <- ifelse(dat=="D", 4, dat)
  dat <- ifelse(is.na(dat), NA, dat)
  dat <- dat %>% factor(., levels = c("1","2","3","4"))
  dat %>% return()
}
UTIdata.df[,56] %<>% gene17.f(.) 

## Select data
UTIdata.df %<>% dplyr::filter(.,!is.na(Hospital_day)) %>% 
  dplyr::filter(.,!is.na(Recurrent)) %>% 
  dplyr::filter(.,!(Hospital_day==0)) %>% 
  dplyr::filter(.,!(ER_H_OPD=="OPD"))

## remove variables with mising rate higher than 0.4 
rm_gene_name.v <- c() 
for (name in paste0("Gene",c(1:16,18:33))) {
  na.prop <- UTIdata.df[,name] %>% table(.,useNA="always") %>% prop.table() %>% as.matrix() %>% .[3]
  if(!is.na(na.prop)){
    if(na.prop>0.4){
      rm_gene_name.v <- c(rm_gene_name.v, name)
    }
  }
}
rm_gene_name.v <- c("Gene1",rm_gene_name.v)

## remove anti variables with mising rate higher than 0.4  
rm_anti_name.v <- c()
for (name in paste0("Anti",c(1:37))) {
  na.prop <- UTIdata.df[,name] %>% table(.,useNA="always") %>% prop.table() %>% as.matrix() %>% .[4]
  if(!is.na(na.prop)){
    if(na.prop>0.4){
      rm_anti_name.v <- c(rm_anti_name.v, name)
    }
  }
}
rm_anti_name.v <- c(rm_anti_name.v,"Anti6","Anti15","Anti17","Anti20","Anti21","Anti22","Anti24","Anti26","Anti27","Anti28",
                    "Anti29","Anti30", "Anti31","Anti32","Anti33","Anti34","Anti35","Anti36","Anti37")
rm_anti_name.v <- unique(rm_anti_name.v)

rm.vector.name.v <- c("Ecoli_No","Date","hospital_long",
                      "Clinical_features","Symptomatic","UWBC","URBC" )
rm_vector_name.v=which(colnames(UTIdata.df) %in% rm.vector.name.v)
rm_anti_idx.v=which(colnames(UTIdata.df) %in% rm_anti_name.v)
rm_gene_idx.v=which(colnames(UTIdata.df) %in% rm_gene_name.v)

with_NA_data <- UTIdata.df[,-c(rm_vector_name.v, rm_anti_idx.v, rm_gene_idx.v)] 
with_NA_data$ER_H_OPD <- ifelse(with_NA_data$ER_H_OPD=="ER-H",1,0)
with_NA_data_R <- with_NA_data %>% dplyr::filter(.,Recurrent==1)
with_NA_data_NR <- with_NA_data %>% dplyr::filter(.,Recurrent==0)
# write.csv(with_NA_data,"C:/Users/88692/Desktop/Scientific_Reports_MLPD_RUTI/Stage2_Before_Imputing_RUTI_Data_Without_PatientInfo.csv")

####################################################################
#####   Imputation for missing values                          #####
####################################################################
without_NA_data_NR=bnstruct::knn.impute(as.matrix(with_NA_data_NR))
without_NA_data_R=bnstruct::knn.impute(as.matrix(with_NA_data_R))
all_data = rbind(without_NA_data_NR, without_NA_data_R) %>% as.data.frame()
recurrent <- c(rep(0, nrow(without_NA_data_NR)), rep(1, nrow(without_NA_data_R)))
all_data <- data.frame(all_data, Recurrent_1=recurrent)
# write.csv(all_data,"C:/Users/88692/Desktop/Scientific_Reports_MLPD_RUTI/Stage2_Imputed_RUTI_Data_Without_PatientInfo.csv")

#################################################
#####       Modeling                        ##### 
#################################################
all_data=read.csv("C:/Users/88692/Desktop/Scientific_Reports_MLPD_RUTI/Stage2_Imputed_RUTI_Data_Without_PatientInfo.csv")
all_data=all_data[,-1]
all_data <- all_data[,-which(colnames(all_data) %in% "Recurrent")]

anti.combine.idx <- c("Anti4","Anti7","Anti18","Anti25")
anti.remove.idx <- c("Anti1","Anti19")
for (name in anti.combine.idx) {
  all_data[,name] <- ifelse(all_data[,name]==2, 1, all_data[,name])
}
all_data <- all_data[,-which(colnames(all_data) %in% anti.remove.idx)]

numeric.idx <- which(colnames(all_data) %in% c("Age","Hospital_day","Pre_hos_2y","Pre_UTI_ER_2y","Pre_UTI_hos_2y","BloodWBC","Creatinine"))
for (i in (1:ncol(all_data))[-numeric.idx]) {
  all_data[,i] <- factor(all_data[,i])
}

all_data[,"UWBC_level"] %<>% ordered(., levels = c(1,2,3,4)) 
all_data[,"URBC_level"] %<>% ordered(., levels = c(0,1,2,3,4)) 
all_data[,"UBact"] %<>% ordered(., levels = c(0,1,2,3,4))
for (j in 2:12) { # for ANTI
  all_data[,j] %<>% ordered(., levels = c(1,2,3)) 
}

## data for Logistic
all_data_glm <- all_data[,-which(colnames(all_data) %in% c("UTI_symptom3",
                                                           "Anti5","Anti8","Anti11","Anti14","Anti23",
                                                           "Dis3","Dis7","Dis8","Dis10","Dis11",
                                                           "Gene4",
                                                           "URBC_level","UTI_symptom3","UTI_symptom8","UTI_symptom6"))]
stage2_RUTI.list <- all_data_glm %>% dplyr::filter(.,Recurrent_1==1)%>% .[1:100,] %>% split(., 1:5)
stage2_UTI.list <- all_data_glm %>% dplyr::filter(.,Recurrent_1==0)  %>% split(., 1:5)
res.idx <- ncol(all_data_glm)

##########################################################
#####             Logistic for Recurrent             #####
##########################################################
glm.f <- function(thres){
  fold_train_accuracy.v=c()
  fold_train_sensitivity.v=c()
  fold_train_specificity.v=c()
  fold_test_accuracy.v=c()
  fold_test_sensitivity.v=c()
  fold_test_specificity.v=c()
  for (t in 1:5) {
    
    testdata = rbind(stage2_RUTI.list[[t]], stage2_UTI.list[[t]])
    testy = testdata[,"Recurrent_1"]
    testy %>% table
    
    train.RUTI.data = do.call(rbind,stage2_RUTI.list[(1:5)[-t]])
    set.seed(1)
    train.UTI.data = do.call(rbind,stage2_UTI.list[(1:5)[-t]]) %>% split(., 1:7)
    
    train_accuracy.v=c()
    train_sensitivity.v=c()
    train_specificity.v=c()
    test.pre.df = data.frame()
    for (k in 1:7) {
      
      traindata <- rbind(train.RUTI.data,train.UTI.data[[k]])
      apply(traindata, 2, table)
      trainy <- traindata$Recurrent_1
      set.seed(1)
      glm.uti=glm(trainy~., data=traindata[,-c(res.idx)], family=binomial)
      
      set.seed(1)
      glm_yhat <- predict(glm.uti, type = "response")
      glm_pred <- ifelse(glm_yhat>0.5,1,0) %>% as.numeric()
      train_accuracy.v <- mean(glm_pred == trainy) %>% c(train_accuracy.v, .)
      result.m <- as.matrix(table(glm_pred, trainy)) 
      train_sensitivity.v <- (result.m[2,2]/(result.m[2,2]+result.m[1,2])) %>% c(train_sensitivity.v, .)
      train_specificity.v <- (result.m[1,1]/(result.m[1,1]+result.m[2,1])) %>% c(train_specificity.v, .)
      
      set.seed(1)
      glm_yhat = predict(glm.uti, newdata=testdata[,-c(res.idx)], type = "response")
      glm_pred <- ifelse(glm_yhat>0.5,1,0) %>% as.numeric()
      
      if(k==1){
        test.pre.df = rep(0, nrow(testdata))
      }
      test.pre.df = cbind(test.pre.df, as.numeric(as.matrix(glm_pred)))
    }
    test.pre.df = test.pre.df[,-1]
    
    maxcount.f = function(dat){
      if(sum(dat)>=thres){
        result=1
      }else{
        result=0
      }
      return(result)
    }
    result.pre = apply(test.pre.df, 1, maxcount.f)
    fold_test_accuracy.v <- mean(result.pre == testy) %>% c(fold_test_accuracy.v, .)
    result.m <- as.matrix(table(result.pre, testy)) 
    fold_test_sensitivity.v <- (result.m[2,2]/(result.m[2,2]+result.m[1,2])) %>% c(fold_test_sensitivity.v, .)
    fold_test_specificity.v <- (result.m[1,1]/(result.m[1,1]+result.m[2,1])) %>% c(fold_test_specificity.v, .)
    fold_train_accuracy.v %<>% c(., mean(train_accuracy.v))
    fold_train_sensitivity.v %<>% c(., mean(train_sensitivity.v))
    fold_train_specificity.v %<>% c(., mean(train_specificity.v))
  }
  glm_result.v=c(round(mean(fold_train_accuracy.v),3), round(sd(fold_train_accuracy.v),3), 
                 round(mean(fold_train_sensitivity.v),3), round(sd(fold_train_sensitivity.v),3), 
                 round(mean(fold_train_specificity.v),3), round(sd(fold_train_specificity.v),3),
                 round(mean(fold_test_accuracy.v),3), round(sd(fold_test_accuracy.v),3), 
                 round(mean(fold_test_sensitivity.v),3), round(sd(fold_test_sensitivity.v),3), 
                 round(mean(fold_test_specificity.v),3), round(sd(fold_test_specificity.v),3))
  return(glm_result.v)
}
glm.result.v <- glm.f(thres = 4)

## data for RandomForest, Tree
stage2_RUTI.list <- all_data %>% dplyr::filter(.,Recurrent_1==1)%>% .[1:100,] %>% split(., 1:5)
stage2_UTI.list <- all_data %>% dplyr::filter(.,Recurrent_1==0)  %>% split(., 1:5)
res.idx <- ncol(all_data)
##############################################################
#####             randomForest for Recurrent             #####
##############################################################
set.seed(1)
rf_model=randomForest(Recurrent_1~., data=all_data, importance=TRUE, proximity=TRUE)

# Figure 3. 
# Variable importance plot of the second stage RF analysis in percentage of mean decrease accuracy for the factors. 
varImpPlot(rf_model)
rf.f <- function(nodesize.n,thres){
  fold_train_accuracy.v=c()
  fold_train_sensitivity.v=c()
  fold_train_specificity.v=c()
  fold_test_accuracy.v=c()
  fold_test_sensitivity.v=c()
  fold_test_specificity.v=c()
  for (t in 1:5) {
    
    testdata = rbind(stage2_RUTI.list[[t]], stage2_UTI.list[[t]])
    testy = testdata[,"Recurrent_1"] %>% factor()
    
    train.RUTI.data = do.call(rbind,stage2_RUTI.list[(1:5)[-t]])
    set.seed(1)
    train.UTI.data = do.call(rbind,stage2_UTI.list[(1:5)[-t]]) %>% split(., 1:7)
    
    train_accuracy.v=c()
    train_sensitivity.v=c()
    train_specificity.v=c()
    test.pre.df = data.frame()
    for (k in 1:7) {
      
      traindata <- rbind(train.RUTI.data,train.UTI.data[[k]])
      trainy <- traindata$Recurrent_1 %>% factor()
      set.seed(1)
      rf_model=randomForest(trainy~., data=traindata[,-c(res.idx)], 
                          importance=TRUE, proximity=TRUE, nodesize=nodesize.n)
      #varImpPlot(rf_model)
      set.seed(1)
      rf_yhat <- predict(rf_model, newdata=traindata[,-c(res.idx)])
      train_accuracy.v <- mean(rf_yhat == trainy) %>% c(train_accuracy.v, .)
      result.m <- as.matrix(table(rf_yhat, trainy)) 
      train_sensitivity.v <- (result.m[2,2]/(result.m[2,2]+result.m[1,2])) %>% c(train_sensitivity.v, .)
      train_specificity.v <- (result.m[1,1]/(result.m[1,1]+result.m[2,1])) %>% c(train_specificity.v, .)
      
      set.seed(1)
      rf_yhat = predict(rf_model, newdata=testdata[,-c(res.idx)])
      if(k==1){
        test.pre.df = rep(0, nrow(testdata))
      }
      test.pre.df = cbind(test.pre.df, as.numeric(as.matrix(rf_yhat)))
    }
    test.pre.df = test.pre.df[,-1]
    
    maxcount.f = function(dat){
      if(sum(dat)>=thres){
        result=1
      }else{
        result=0
      }
      return(result)
    }
    result.pre = apply(test.pre.df, 1, maxcount.f)
    fold_test_accuracy.v <- mean(result.pre == testy) %>% c(fold_test_accuracy.v, .)
    result.m <- as.matrix(table(result.pre, testy)) 
    fold_test_sensitivity.v <- (result.m[2,2]/(result.m[2,2]+result.m[1,2])) %>% c(fold_test_sensitivity.v, .)
    fold_test_specificity.v <- (result.m[1,1]/(result.m[1,1]+result.m[2,1])) %>% c(fold_test_specificity.v, .)
    fold_train_accuracy.v %<>% c(., mean(train_accuracy.v))
    fold_train_sensitivity.v %<>% c(., mean(train_sensitivity.v))
    fold_train_specificity.v %<>% c(., mean(train_specificity.v))
  }
  rf_result.v=c(round(mean(fold_train_accuracy.v),3), round(sd(fold_train_accuracy.v),3), 
                round(mean(fold_train_sensitivity.v),3), round(sd(fold_train_sensitivity.v),3), 
                round(mean(fold_train_specificity.v),3), round(sd(fold_train_specificity.v),3),
                round(mean(fold_test_accuracy.v),3), round(sd(fold_test_accuracy.v),3), 
                round(mean(fold_test_sensitivity.v),3), round(sd(fold_test_sensitivity.v),3), 
                round(mean(fold_test_specificity.v),3), round(sd(fold_test_specificity.v),3))
  return(rf_result.v)
}
rf.result.v <- rf.f(nodesize.n = 5, thres = 5)

######################################################
#####             TREE for Recurrent             #####
######################################################
set.seed(1)
tree_model=rpart(Recurrent_1~., data=all_data) 

# Figure 4. 
# The decision rules of the DT analysis for development of RUTI after hospitalization.
rpart.plot::rpart.plot(tree_model, extra=1, digits=2, split.col="red", cex=1)
yhat.tree <- predict(tree_model, newdata=all_data, type = c("class"))
accuracy <- mean(yhat.tree == all_data$Recurrent_1)
result.m <- as.matrix(table(yhat.tree, all_data$Recurrent_1)) 
sensitivity <- (result.m[2,2]/(result.m[2,2]+result.m[1,2])) 
specificity <- (result.m[1,1]/(result.m[1,1]+result.m[2,1]))

tree.f <- function(minsplit.n,thres){
  fold_train_accuracy.v=c()
  fold_train_sensitivity.v=c()
  fold_train_specificity.v=c()
  fold_test_accuracy.v=c()
  fold_test_sensitivity.v=c()
  fold_test_specificity.v=c()
  for (t in 1:5) {
    
    testdata = rbind(stage2_RUTI.list[[t]], stage2_UTI.list[[t]])
    testy = testdata[,"Recurrent_1"] %>% factor()
    
    train.RUTI.data = do.call(rbind,stage2_RUTI.list[(1:5)[-t]])
    set.seed(1)
    train.UTI.data = do.call(rbind,stage2_UTI.list[(1:5)[-t]]) %>% split(., 1:7)
    
    train_accuracy.v=c()
    train_sensitivity.v=c()
    train_specificity.v=c()
    test.pre.df = data.frame()
    for (k in 1:7) {
      
      traindata <- rbind(train.RUTI.data,train.UTI.data[[k]])
      trainy <- traindata$Recurrent_1 %>% factor()
      
      set.seed(1)
      tree_model=rpart(trainy~., data=traindata[,-c(res.idx)], minsplit = minsplit.n) 
      set.seed(1)
      tree_yhat = predict(tree_model, newdata=traindata[,-c(res.idx)], type = c("class"))
      
      train_accuracy.v <- mean(tree_yhat == trainy) %>% c(train_accuracy.v, .)
      result.m <- as.matrix(table(tree_yhat, trainy)) 
      train_sensitivity.v <- (result.m[2,2]/(result.m[2,2]+result.m[1,2])) %>% c(train_sensitivity.v, .)
      train_specificity.v <- (result.m[1,1]/(result.m[1,1]+result.m[2,1])) %>% c(train_specificity.v, .)
      
      set.seed(20)
      tree_yhat = predict(tree_model, newdata=testdata[,-c(res.idx)], type = c("class"))
      if(k==1){
        test.pre.df = rep(0, nrow(testdata))
      }
      test.pre.df = cbind(test.pre.df, as.numeric(as.matrix(tree_yhat)))
    }
    test.pre.df = test.pre.df[,-1]
    
    maxcount.f = function(dat){
      if(sum(dat)>=thres){
        result=1
      }else{
        result=0
      }
      return(result)
    }
    
    result.pre = apply(test.pre.df, 1, maxcount.f)
    fold_test_accuracy.v <- mean(result.pre == testy) %>% c(fold_test_accuracy.v, .)
    result.m <- as.matrix(table(result.pre, testy)) 
    fold_test_sensitivity.v <- (result.m[2,2]/(result.m[2,2]+result.m[1,2])) %>% c(fold_test_sensitivity.v, .)
    fold_test_specificity.v <- (result.m[1,1]/(result.m[1,1]+result.m[2,1])) %>% c(fold_test_specificity.v, .)
    fold_train_accuracy.v %<>% c(., mean(train_accuracy.v))
    fold_train_sensitivity.v %<>% c(., mean(train_sensitivity.v))
    fold_train_specificity.v %<>% c(., mean(train_specificity.v))
  }
  tree_result.v=c(round(mean(fold_train_accuracy.v),3), round(sd(fold_train_accuracy.v),3), 
                    round(mean(fold_train_sensitivity.v),3), round(sd(fold_train_sensitivity.v),3), 
                    round(mean(fold_train_specificity.v),3), round(sd(fold_train_specificity.v),3),
                    round(mean(fold_test_accuracy.v),3), round(sd(fold_test_accuracy.v),3), 
                    round(mean(fold_test_sensitivity.v),3), round(sd(fold_test_sensitivity.v),3), 
                    round(mean(fold_test_specificity.v),3), round(sd(fold_test_specificity.v),3))
  return(tree_result.v)
}
tree.result.v=tree.f(minsplit.n = 5, thres = 4)

## Table 4. 
## Comparison of the performance in RUTI prediction models after hospitalization for UTI through 5-fold cross validation (sample size = 809). 
all_result.df=rbind(glm.result.v[7:12],tree.result.v[7:12],rf.result.v[7:12])
result_colnames1 <- c(rep(c("Accuracy","Sensitivity","Specificity"),each=2))
result_colnames2 <- c(rep(c("Mean","Standard deviation"),3))
all_result.df <- rbind(result_colnames1, result_colnames2, all_result.df)

###############################################################################
## ?p??table 2?Mtable3?ƭ?
temp_UTIdata <- read.csv("D:/UTIdata/Scientific_Reports_MLPD_RUTI/Stage2_Before_Imputing_RUTI_Data_Without_PatientInfo.csv")
temp_UTIdata <- as.data.frame(temp_UTIdata)

table2_data=temp_UTIdata[,c('Gene17', "Gene2", "Gene3", "Gene4", "Gene5", 
                            "Gene6", "Gene7", "Gene8", "Gene9", "Gene10", 
                            "Gene11", "Gene12", "Gene13", "Gene14", "Gene15", "Gene16",
                            "Anti1","Anti2","Anti4","Anti5","Anti7","Anti8","Anti11",
                            "Anti13","Anti14","Anti18","Anti19","Anti23","Anti25")]
table2_RUTI=temp_UTIdata[,'Recurrent']

# gene 17
colnames(table2_data)[1]
N1=length(na.omit(table2_data[table2_RUTI==0,1]))
N2=length(na.omit(table2_data[table2_RUTI==1,1]))
n1=sum(table(table2_data[table2_RUTI==0,1]))
n2=sum(table(table2_data[table2_RUTI==1,1]))
gene17.m=cbind(rep('gene17',5),
               c('n = ',as.numeric(table(table2_data[table2_RUTI==0,1]))),
               c(N1,round((table(table2_data[table2_RUTI==0,1])/n1)*100)),
               c('n = ',as.numeric(table(table2_data[table2_RUTI==1,1]))),
               c(N2,round((table(table2_data[table2_RUTI==1,1])/n2)*100)),
               round(fisher.test(matrix(cbind(table(table2_data[table2_RUTI==0,1]),
                                              table(table2_data[table2_RUTI==1,1])), nrow = 4))$p.value,4))

# gene
colnames(table2_data)[2:16]
gene.m=matrix(NA,nrow = 15*2,ncol = 6)
for (i in 2:16) {
  idx_1=table2_data[,i]
  N1=length(na.omit(table2_data[table2_RUTI==0,i]))
  N2=length(na.omit(table2_data[table2_RUTI==1,i]))
  n1=length(na.omit(table2_data[idx_1==1 & table2_RUTI==0,i]))
  n2=length(na.omit(table2_data[idx_1==1 & table2_RUTI==1,i]))
  pvalue=fisher.test(matrix(c(n1, N1-n1, n2, N2-n2),nrow=2))$p.value
  gene.m[c(2*(i-1)-1,2*(i-1)),]=rbind(t(matrix(c(colnames(table2_data)[i],
                                                 'n = ', N1,
                                                 'n = ', N2,
                                                 round(pvalue,4)))),
                                      t(matrix(c(colnames(table2_data)[i],
                                                 n1, round((n1/N1)*100),
                                                 n2, round((n2/N2)*100),
                                                 round(pvalue,4)))))
}

## gene2 or gene3
temp.m=cbind(table2_data[,c(2,3)],rep(0,nrow(table2_data)),table2_RUTI)
idx_1=rowSums(temp.m[,1:3],na.rm = T)
idx_2=ifelse(idx_1!=0,1,0)
idx_3=ifelse(!is.na(temp.m[,1]) | !is.na(temp.m[,2]),1,0)
N1=sum(table2_RUTI==0 & idx_3==1)
N2=sum(table2_RUTI==1 & idx_3==1)
n1=sum(idx_2==1 & table2_RUTI==0)
n2=sum(idx_2==1 & table2_RUTI==1)
pvalue=fisher.test(matrix(c(n1, N1-n1, n2, N2-n2),nrow=2))$p.value


# anti
colnames(table2_data)[17:29]
anti.m=matrix(NA,nrow = 13*4,ncol = 6)
for (i in 17:29) {
  idx_1=table2_data[,i]
  N1=length(na.omit(table2_data[table2_RUTI==0,i]))
  N2=length(na.omit(table2_data[table2_RUTI==1,i]))
  n10=length(na.omit(table2_data[idx_1==1 & table2_RUTI==0,i])) 
  n20=length(na.omit(table2_data[idx_1==2 & table2_RUTI==0,i])) 
  n30=length(na.omit(table2_data[idx_1==3 & table2_RUTI==0,i])) 
  n11=length(na.omit(table2_data[idx_1==1 & table2_RUTI==1,i])) 
  n21=length(na.omit(table2_data[idx_1==2 & table2_RUTI==1,i])) 
  n31=length(na.omit(table2_data[idx_1==3 & table2_RUTI==1,i])) 
  pvalue=round(chisq.test(matrix(c(n10, n20, n30, n11, n21, n31), nrow = 3))$p.value,4)
  anti.m[(4*(i-17)+1),1]=colnames(table2_data)[i]
  anti.m[(4*(i-17)+1),c(2,4)]='n = '
  anti.m[(4*(i-17)+1),3]=N1
  anti.m[(4*(i-17)+1),5]=N2
  anti.m[(4*(i-17)+1),6]=pvalue
  anti.m[(4*(i-17)+2):(4*(i-17)+4),2:5]=cbind(c(n10,n20,n30),
                                              c(round((n10/N1)*100),round((n20/N1)*100),round((n30/N1)*100)),
                                              c(n11,n21,n31),
                                              c(round((n11/N2)*100),round((n21/N2)*100),round((n31/N2)*100)))
  
}

## Table 2. 
# Bacterial characteristics related to UTI and RUTI (sample size = 809) used in the second stage analysis.
table2_result.m=rbind(gene17.m,gene.m)
s1=paste0(table2_result.m[,2],' (',table2_result.m[,3],')')
s2=paste0(table2_result.m[,4],' (',table2_result.m[,5],')')
table2_result.m=cbind(table2_result.m[,1],s1,s2,table2_result.m[,6],table2_result.m[,c(2:5)])

## Table 3.
## Antimicrobial susceptibility of bacterial pathogens related to UTI and RUTI (sample size = 809) used in the second stage analysis.
table3_result.m=anti.m
s1=paste0(table3_result.m[,2],' (',table3_result.m[,3],')')
s2=paste0(table3_result.m[,4],' (',table3_result.m[,5],')')
table3_result.m=cbind(table3_result.m[,1],s1,s2,table3_result.m[,6],table3_result.m[,c(2:5)])

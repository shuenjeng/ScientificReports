# The code for the paper
# Machine learning to predict the development of recurrent urinary tract infection related to single uropathogen, Escherichia coli
#
#
#Shuen-Lin Jeng1, Zi-Jing Huang2, Deng-Chi Yang3, Ching-Hao Teng4*, Ming-Cheng Wang5*
#
#1Department of Statistics, Institute of Data Science, and Center for Innovative FinTech Business Models, National Cheng Kung University, Tainan, Taiwan
#2Department of Statistics, National Cheng Kung University, Tainan, Taiwan
#3Department of Geriatrics and Gerontology, National Cheng Kung University Hospital, College of Medicine, National Cheng Kung University, Tainan, Taiwan.
#4Institute of Molecular Medicine, College of Medicine, National Cheng Kung University, Tainan, Taiwan; Institute of Basic Medical Sciences, College of Medicine, National Cheng Kung University, Tainan, Taiwan; Center of Infectious Disease and Signaling Research, National Cheng Kung University, Tainan, Taiwan.
#5Division of Nephrology, Department of Internal Medicine, National Cheng Kung University Hospital, College of Medicine, National Cheng Kung University, Tainan, Taiwan.
#
#================================================================================================================================================
library(readxl)
library(magrittr)
library(randomForest)
library(rpart)
library(rpart.plot)
options(stringsAsFactors = FALSE)

#===========================================================
# First stage Analysis

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

UTIdata.df$ER_H_OPD %<>% ifelse(. %in% c("ER-NO H","ER-H (DEAD)","ER-H"), "ER", .) %>% factor()
UTIdata.df$Bacterial_Name %<>% ifelse(. %in% c("Klebsiella pneumoniae"), NA, .) %>% as.factor(.)
UTIdata.df$Recurrent %<>% as.factor(.)
UTIdata.df$Dead %<>% ifelse(.==9, NA, .) %>% as.factor()
UTIdata.df$Gender %<>% as.factor(.)
UTIdata.df$UTI_pos %<>% ifelse(. %in% c("X","x","0"), NA, .) %>% as.factor()
UTIdata.df$Pre_hos_2y %<>% ifelse(. %in% c("x"), NA, .) %>% as.numeric()
UTIdata.df$Pre_UTI_ER_2y %<>% ifelse(. %in% c("x"), NA, .) %>% as.numeric()
UTIdata.df$Pre_UTI_hos_2y %<>% ifelse(. %in% c("x"), NA, .) %>% as.numeric()
UTIdata.df$BloodWBC %<>% ifelse(. %in% c("x"), NA, .) %>% as.numeric()
UTIdata.df$Creatinine %<>% ifelse(. %in% c("x"), NA, .) %>% as.numeric()
UTIdata.df$UBact %<>% ifelse(. %in% c("x"), NA, .) %>% as.factor()
UTIdata.df$UBact %<>% ordered(., levels = c(0,1,2,3,4))
UTIdata.df$Age %<>% ifelse(.==767,NA,.)
UTIdata.df[,'four_disease_group'] = with(UTIdata.df,ifelse(Dis5==1|Dis6==1|Dis7==1|Dis9==1,1,0)) # ?c???D?????\?????`
UTIdata.df[,'four_disease_group'] %<>% factor(.)

UTIdata.df$UTI_symptom1 %<>% ifelse(. %in% c(2,11), NA, .) 
UTIdata.df$UTI_symptom5 %<>% ifelse(. %in% c("x"), NA, .)
UTIdata.df$UTI_symptom8 %<>% ifelse(. == 19.9, NA, .)
UTIdata.df[,c(97:104)] <- apply(UTIdata.df[,c(97:104)], 2, as.factor)

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
UTIdata.df$UWBC_level <- sapply(1:length(UTIdata.df$UWBC), WBC_level.f, dat=UTIdata.df$UWBC)
UTIdata.df$UWBC_level %<>%  ordered(., levels = c(1,2,3,4))

UTIdata.df$URBC %<>% ifelse(. %in% c("x"), NA, .)
UTIdata.df$URBC %<>% ifelse(. %in% c("numerous"), 1001, .)
UTIdata.df$URBC %<>% ifelse(. %in% c("100+"), 101, .)
UTIdata.df$URBC %<>% ifelse(. %in% c("1000+"), 1001, .)
UTIdata.df$URBC %<>% as.numeric(.)
UTIdata.df$URBC_level <- sapply(1:length(UTIdata.df$URBC), RBC_level.f, dat=UTIdata.df$URBC)
UTIdata.df$URBC_level %<>% ordered(., levels = c(0,1,2,3,4))

#############################################################################
#####                Choose variable for fitting model                  #####
#####  Select Non-Hospital Patient and Convert to categorical variable  #####
#############################################################################
later_var.v <- c("Gender","Age",paste0("Dis",c(1:10,12)),
                 "Pre_hos_2y","Pre_UTI_ER_2y","Pre_UTI_hos_2y",
                 paste0("UTI_symptom",c(1:2,4:8)), 
                 "BloodWBC","Creatinine","UWBC_level","URBC_level",
                 "UBact","ER_H_OPD",'four_disease_group')
analysis_data.df <- data.frame(UTIdata.df[,later_var.v], 
                               Recurrent=UTIdata.df$Recurrent) %>% na.omit()
## Change to categorical data
for (name in c("Gender","four_disease_group",paste0("Dis",c(1:10,12)),paste0("UTI_symptom",c(1:2,4:8)))) {
  analysis_data.df[,name] <- factor(analysis_data.df[,name])
}

## Select data
analysis_data.df=analysis_data.df[analysis_data.df$ER_H_OPD%in%c('ER','OPD'),]
analysis_data.df$ER_H_OPD %<>% factor(.)

## 5-fold split
stage1_RUTI.list <- analysis_data.df %>% dplyr::filter(., Recurrent==1) %>% split(., 1:5)
stage1_UTI.list <- analysis_data.df %>% dplyr::filter(., Recurrent==0) %>% split(., 1:5)

##########################################################
#####             Logistic for Recurrent             #####
##########################################################
glm.f <- function(threshold,RUTI.list,UTI.list,split_num){ 
  fold_train_accuracy.v=c()
  fold_train_sensitivity.v=c()
  fold_train_specificity.v=c()
  fold_test_accuracy.v=c()
  fold_test_sensitivity.v=c()
  fold_test_specificity.v=c()
  for (t in 1:5) {
    
    testdata = rbind(RUTI.list[[t]], UTI.list[[t]])
    testy = testdata[,"Recurrent"]
    
    train_RUTI_data.df = do.call(rbind, RUTI.list[(1:5)[-t]])
    set.seed(1)
    train_UTI_data.df = do.call(rbind, UTI.list[(1:5)[-t]]) %>% split(., 1:split_num)
    
    train_accuracy.v=c()
    train_sensitivity.v=c()
    train_specificity.v=c()
    test.pre.df = data.frame()
    for (k in 1:split_num) {
      
      traindata <- rbind.data.frame(train_RUTI_data.df,train_UTI_data.df[[k]])
      trainy <- traindata$Recurrent
      set.seed(1)
      res.idx=which(colnames(traindata)=='Recurrent')
      glm_model=glm(trainy~., data=traindata[,-c(res.idx)], family=binomial)
      
      set.seed(1)
      glm_yhat <- predict(glm_model, type = "response")
      glm_pred <- ifelse(glm_yhat>0.5,1,0) %>% as.numeric()
      train_accuracy.v <- mean(glm_pred == trainy) %>% c(train_accuracy.v, .)
      result.m <- as.matrix(table(glm_pred, trainy)) 
      train_sensitivity.v <- (result.m[2,2]/(result.m[2,2]+result.m[1,2])) %>% c(train_sensitivity.v, .)
      train_specificity.v <- (result.m[1,1]/(result.m[1,1]+result.m[2,1])) %>% c(train_specificity.v, .)
      
      set.seed(1)
      res.idx=which(colnames(traindata)=='Recurrent')
      glm_yhat = predict(glm_model, newdata=testdata[,-c(res.idx)], type = "response")
      glm_pred <- ifelse(glm_yhat>0.5,1,0) %>% as.numeric()
      
      if(k==1){
        test.pre.df = rep(0, nrow(testdata))
      }
      test.pre.df = cbind(test.pre.df, as.numeric(as.matrix(glm_pred)))
    }
    test.pre.df = test.pre.df[,-1]
    
    maxcount.f = function(dat){
      if(sum(dat)>=threshold){
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
glm.stage1_RUTI.list <- analysis_data.df[,-5] %>% dplyr::filter(.,Recurrent==1) %>% split(., 1:5)
glm.stage1_UTI.list <- analysis_data.df[,-5] %>% dplyr::filter(.,Recurrent==0) %>% split(., 1:5)
glm.result.v=glm.f(threshold=2,
                   RUTI.list=glm.stage1_RUTI.list, 
                   UTI.list=glm.stage1_UTI.list,
                   split_num=5)

##############################################################
#####             randomForest for Recurrent             #####
##############################################################
set.seed(1)
rf_model=randomForest(Recurrent~., 
                      data=analysis_data.df, 
                      importance=TRUE, proximity=TRUE)

## Figure 1. 
## Variable importance plot of the first stage RF analysis in percentage of mean decrease accuracy for the factors.
varImpPlot(rf_model,n.var=30)
rf_yhat <- predict(rf_model, newdata=analysis_data.df)
accuracy <- mean(rf_yhat == analysis_data.df$Recurrent)
result.m <- as.matrix(table(rf_yhat, analysis_data.df$Recurrent)) 
sensitivity <- (result.m[2,2]/(result.m[2,2]+result.m[1,2])) 
specificity <- (result.m[1,1]/(result.m[1,1]+result.m[2,1]))
rf.f <- function(nodesize.n,threshold,RUTI.list,UTI.list,split_num=5,ntree.n=500){
  fold_train_accuracy.v=c()
  fold_train_sensitivity.v=c()
  fold_train_specificity.v=c()
  fold_test_accuracy.v=c()
  fold_test_sensitivity.v=c()
  fold_test_specificity.v=c()
  for (t in 1:5) {
    
    testdata = rbind(RUTI.list[[t]], UTI.list[[t]])
    testy = testdata[,"Recurrent"]
    
    train_RUTI_data.df = do.call(rbind,RUTI.list[(1:5)[-t]])
    set.seed(1)
    train_UTI_data.df = do.call(rbind,UTI.list[(1:5)[-t]]) %>% split(., 1:split_num)
    
    train_accuracy.v=c()
    train_sensitivity.v=c()
    train_specificity.v=c()
    test.pre.df = data.frame()
    for (k in 1:split_num) {
      
      traindata <- rbind(train_RUTI_data.df,train_UTI_data.df[[k]])
      trainy <- traindata$Recurrent
      set.seed(1)
      res.idx=which(colnames(traindata)=='Recurrent')
      rf_model=randomForest(trainy~., data=traindata[,-c(res.idx)], 
                            importance=TRUE, proximity=TRUE, 
                            nodesize=nodesize.n,ntree=ntree.n)
      set.seed(1)
      res.idx=which(colnames(traindata)=='Recurrent')
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
      if(sum(dat)>=threshold){
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
rf.result.v=rf.f(nodesize.n=20,
                 threshold=3,
                 RUTI.list=stage1_RUTI.list,
                 UTI.list=stage1_UTI.list,
                 split_num=5,
                 ntree.n=500)


######################################################
#####             TREE for Recurrent             #####
######################################################
set.seed(1)
tree_model=rpart(Recurrent~., data=analysis_data.df)

## Figure 2.
## The decision rules of the DT analysis for development of RUTI in the clinical visit. (sample size = 963). 
rpart.plot::rpart.plot(tree_model, extra=1, digits=2, split.col="red", cex=1)
tree_yhat <- predict(tree_model, newdata=analysis_data.df, type = c("class"))
accuracy <- mean(tree_yhat == analysis_data.df$Recurrent)
result.m <- as.matrix(table(tree_yhat, analysis_data.df$Recurrent)) 
sensitivity <- (result.m[2,2]/(result.m[2,2]+result.m[1,2])) 
specificity <- (result.m[1,1]/(result.m[1,1]+result.m[2,1]))
tree.f <- function(minsplit.n,threshold,RUTI.list,UTI.list,split_num){
  fold_train_accuracy.v=c()
  fold_train_sensitivity.v=c()
  fold_train_specificity.v=c()
  fold_test_accuracy.v=c()
  fold_test_sensitivity.v=c()
  fold_test_specificity.v=c()
  for (t in 1:5) {
    
    testdata = rbind(RUTI.list[[t]], UTI.list[[t]])
    testy = testdata[,"Recurrent"]
    
    train_RUTI_data.df = do.call(rbind,RUTI.list[(1:5)[-t]])
    set.seed(20)
    train_UTI_data.df = do.call(rbind,UTI.list[(1:5)[-t]]) %>% split(., 1:split_num)
    
    train_accuracy.v=c()
    train_sensitivity.v=c()
    train_specificity.v=c()
    test.pre.df = data.frame()
    for (k in 1:split_num) {
      
      traindata <- rbind(train_RUTI_data.df,train_UTI_data.df[[k]])
      trainy <- traindata$Recurrent
      
      set.seed(20)
      res.idx=which(colnames(traindata)=='Recurrent')
      tree_model=rpart(trainy~., data=traindata[,-c(res.idx)], minsplit = minsplit.n) 
      # rpart.plot::rpart.plot(tree_model, extra=1, digits=2, split.col="red", cex=0.7)
      set.seed(20)
      res.idx=which(colnames(traindata)=='Recurrent')
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
      if(sum(dat)>=threshold){
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
tree.result.v=tree.f(minsplit.n=5,
                     threshold=3,
                     RUTI.list=stage1_RUTI.list,
                     UTI.list=stage1_UTI.list,
                     split_num=5)

## Table 3. 
## Comparison of the performance in RUTI prediction models in the clinical visit through 5-fold cross validation (sample size = 963).
all_result.df=rbind(glm.result.v[7:12],tree.result.v[7:12],rf.result.v[7:12])
result_colnames1 <- c(rep(c("Accuracy","Sensitivity","Specificity"),each=2))
result_colnames2 <- c(rep(c("Mean","Standard deviation"),3))
all_result.df <- rbind(result_colnames1, result_colnames2, all_result.df)

###############################################################################
## ???J????data3.xlsx->?p??table 1????
temp_UTIdata <- read_xlsx("D:/UTIdata/Scientific_Reports_MLPD_RUTI_20220707/data3.xlsx",sheet=1)
temp_UTIdata <- as.data.frame(temp_UTIdata)
table1_data=temp_UTIdata[,c(7,6,36,8:19,20,21,22,23:30,32,31,34,33,35)]
colnames(table1_data)[c(8,9,10,12)] # dis5679: Foley, Obstructive.uropathy, Urolithiasis, Neurogenic.bladder
table1_data[,'four_disease_group']=ifelse(rowSums(table1_data[,c(8,9,10,12)])==0,0,1)
colnames(table1_data)[c(19:26)]
for (i in 19:26) {
  table1_data[,i]=as.numeric(table1_data[,i])
}
table1_data[,'any_utisym']=ifelse(rowSums(table1_data[,c(19:26)])==0,0,1)

colnames(temp_UTIdata)[5]
table1_RUTI=temp_UTIdata[,5]

# age
age.m=t(matrix(c('Age (year)',
                 median(table1_data$Age[table1_RUTI==0]),NA,
                 median(table1_data$Age[table1_RUTI==1]),NA,
                 round(wilcox.test(as.numeric(table1_data$Age[table1_RUTI==0]),
                                   as.numeric(table1_data$Age[table1_RUTI==1]))$p.value,4))))
quantile(table1_data$Age[table1_RUTI==0],probs = c(0.25,0.75))
quantile(table1_data$Age[table1_RUTI==1],probs = c(0.25,0.75))



# gender
colnames(table1_data)[2]
for (i in 2) {
  idx_1=table1_data[,i]
  n1=length(table1_data[idx_1==1 & table1_RUTI==0,i])
  n2=length(table1_data[idx_1==1 & table1_RUTI==1,i])
  pvalue=fisher.test(matrix(c(n1, 826-n1, n2, 137-n2),nrow=2))$p.value
  gender.m=t(matrix(c('Gender (male)',
                      n1, round((n1/826)*100),
                      n2, round((n2/137)*100),
                      round(pvalue,4))))
}

# Place_of_collection
colnames(table1_data)[3]
for (i in 3) {
  idx_1=table1_data[,i]
  n1=length(table1_data[idx_1%in%c("ER-H","ER-NO H") & table1_RUTI==0,i])
  n2=length(table1_data[idx_1%in%c("ER-H","ER-NO H") & table1_RUTI==1,i])
  pvalue=fisher.test(matrix(c(n1, 826-n1, n2, 137-n2),nrow=2))$p.value
  PC.m=t(matrix(c('Place of urine sample collection (ED)  (Place_of_collection)',
                  n1, round((n1/826)*100),
                  n2, round((n2/137)*100),
                  round(pvalue,4))))
}

# disease
colnames(table1_data)[4:15]
dis.m=matrix(NA,nrow = 12,ncol = 6)
for (i in 4:15) {
  idx_1=table1_data[,i]
  n1=length(table1_data[idx_1==1 & table1_RUTI==0,i])
  n2=length(table1_data[idx_1==1 & table1_RUTI==1,i])
  pvalue=fisher.test(matrix(c(n1, 826-n1, n2, 137-n2),nrow=2))$p.value
  dis.m[i-3,]=t(matrix(c(colnames(table1_data)[i],
                         n1, round((n1/826)*100),
                         n2, round((n2/137)*100),
                         round(pvalue,4))))
}

# four_disease_group
colnames(table1_data)[32]
for (i in 32) {
  idx_1=table1_data[,i]
  n1=length(table1_data[idx_1==1 & table1_RUTI==0,i])
  n2=length(table1_data[idx_1==1 & table1_RUTI==1,i])
  pvalue=fisher.test(matrix(c(n1, 826-n1, n2, 137-n2),nrow=2))$p.value
  four_dis.m=t(matrix(c(colnames(table1_data)[i],
                        n1, round((n1/826)*100),
                        n2, round((n2/137)*100),
                        round(pvalue,4))))
}

# previous record
colnames(table1_data)[16:18]
pre_hos.m=t(matrix(c('Frequency of hospitalization within 2 years (Pre_hos_2y)',
                     median(as.numeric(table1_data[table1_RUTI==0,16]),na.rm = T), NA,
                     median(as.numeric(table1_data[table1_RUTI==1,16]),na.rm = T), NA,
                     round(wilcox.test(as.numeric(table1_data[table1_RUTI==0,16]),
                                       as.numeric(table1_data[table1_RUTI==1,16]))$p.value,4))))
quantile(na.omit(as.numeric(table1_data[table1_RUTI==0,16])),probs = c(0.25,0.75))
quantile(na.omit(as.numeric(table1_data[table1_RUTI==1,16])),probs = c(0.25,0.75))

pre_uti_er.m=t(matrix(c('Frequency of ED visit within 2 years (Pre_UTI_ER_2y) ',
                        median(as.numeric(table1_data[table1_RUTI==0,17]),na.rm = T), NA,
                        median(as.numeric(table1_data[table1_RUTI==1,17]),na.rm = T), NA,
                        round(wilcox.test(as.numeric(table1_data[table1_RUTI==0,17]),
                                          as.numeric(table1_data[table1_RUTI==1,17]))$p.value,4))))
quantile(na.omit(as.numeric(table1_data[table1_RUTI==0,17])),probs = c(0.25,0.75))
quantile(na.omit(as.numeric(table1_data[table1_RUTI==1,17])),probs = c(0.25,0.75))

pre_uti_hos.m=t(matrix(c('Frequency of UTI within 2 years (Pre_UTI_hos_2y)',
                         median(as.numeric(table1_data[table1_RUTI==0,18]),na.rm = T), NA,
                         median(as.numeric(table1_data[table1_RUTI==1,18]),na.rm = T), NA,
                         round(wilcox.test(as.numeric(table1_data[table1_RUTI==0,18]),
                                           as.numeric(table1_data[table1_RUTI==1,18]))$p.value,4))))
quantile(na.omit(as.numeric(table1_data[table1_RUTI==0,18])),probs = c(0.25,0.75))
quantile(na.omit(as.numeric(table1_data[table1_RUTI==1,18])),probs = c(0.25,0.75))

# any UTI symptom
colnames(table1_data)[33]
for (i in 33) {
  idx_1=table1_data[,i]
  n1=length(table1_data[idx_1==1 & table1_RUTI==0,i])
  n2=length(table1_data[idx_1==1 & table1_RUTI==1,i])
  pvalue=fisher.test(matrix(c(n1, 826-n1, n2, 137-n2),nrow=2))$p.value
  any_utisym.m=t(matrix(c('Any UTI symptom',
                          n1, round((n1/826)*100),
                          n2, round((n2/137)*100),
                          round(pvalue,4))))
}

# UTI.symptoms
colnames(table1_data)[19:26]
utisym.m=matrix(NA,nrow = 8,ncol = 6)
for (i in 19:26) {
  idx_1=table1_data[,i]
  n1=length(table1_data[idx_1==1 & table1_RUTI==0,i])
  n2=length(table1_data[idx_1==1 & table1_RUTI==1,i])
  pvalue=fisher.test(matrix(c(n1, 826-n1, n2, 137-n2),nrow=2))$p.value
  utisym.m[i-18,]=t(matrix(c(colnames(table1_data)[i],
                             n1, round((n1/826)*100),
                             n2, round((n2/137)*100),
                             round(pvalue,4))))
}

# creatinine
colnames(table1_data)[27]
creatinine.m=t(matrix(c('Serum creatinine (mg/dL)',
                        median(as.numeric(table1_data[table1_RUTI==0,27]),na.rm = T), NA,
                        median(as.numeric(table1_data[table1_RUTI==1,27]),na.rm = T), NA,
                        round(wilcox.test(as.numeric(table1_data[table1_RUTI==0,27]),
                                          as.numeric(table1_data[table1_RUTI==1,27]))$p.value,4))))
quantile(na.omit(as.numeric(table1_data[table1_RUTI==0,27])),probs = c(0.25,0.75))
quantile(na.omit(as.numeric(table1_data[table1_RUTI==1,27])),probs = c(0.25,0.75))

# Peak blood WBC count 
colnames(table1_data)[28]
BWBC.m=t(matrix(c('Peak blood WBC count (109/L) (BloodWBC)',
                  median(as.numeric(table1_data[table1_RUTI==0,28]),na.rm = T), NA,
                  median(as.numeric(table1_data[table1_RUTI==1,28]),na.rm = T), NA,
                  round(wilcox.test(as.numeric(table1_data[table1_RUTI==0,28]),
                                    as.numeric(table1_data[table1_RUTI==1,28]))$p.value,4))))
quantile(na.omit(as.numeric(table1_data[table1_RUTI==0,28])),probs = c(0.25,0.75))
quantile(na.omit(as.numeric(table1_data[table1_RUTI==1,28])),probs = c(0.25,0.75))

# Urinary bacterial count
colnames(table1_data)[29]
UBact.m=t(matrix(c('Urinary bacterial count (0~4) (UBact)',
                   median(as.numeric(table1_data[table1_RUTI==0,29]),na.rm = T),NA,
                   median(as.numeric(table1_data[table1_RUTI==1,29]),na.rm = T),NA,
                   round(wilcox.test(as.numeric(table1_data[table1_RUTI==0,29]),
                                     as.numeric(table1_data[table1_RUTI==1,29]))$p.value,4))))
quantile(na.omit(as.numeric(table1_data[table1_RUTI==0,29])),probs = c(0.25,0.75))
quantile(na.omit(as.numeric(table1_data[table1_RUTI==1,29])),probs = c(0.25,0.75))

# Urinary WBC/HPF
colnames(table1_data)[30]
UWBC.m=t(matrix(c('Urinary WBC/HPF (UWBC_level)',
                  median(as.numeric(table1_data[table1_RUTI==0,30]),na.rm = T),NA,
                  median(as.numeric(table1_data[table1_RUTI==1,30]),na.rm = T),NA,
                  round(wilcox.test(as.numeric(table1_data[table1_RUTI==0,30]),
                                    as.numeric(table1_data[table1_RUTI==1,30]))$p.value,4))))
quantile(na.omit(as.numeric(table1_data[table1_RUTI==0,30])),probs = c(0.25,0.75))
quantile(na.omit(as.numeric(table1_data[table1_RUTI==1,30])),probs = c(0.25,0.75))

# Urinary RBC/HPF 
colnames(table1_data)[31]
URBC.m=t(matrix(c('Urinary RBC/HPF (URBC_level)',
                  median(as.numeric(table1_data[table1_RUTI==0,31]),na.rm = T),NA,
                  median(as.numeric(table1_data[table1_RUTI==1,31]),na.rm = T),NA,
                  round(wilcox.test(as.numeric(table1_data[table1_RUTI==0,31]),
                                    as.numeric(table1_data[table1_RUTI==1,31]))$p.value,4))))
quantile(na.omit(as.numeric(table1_data[table1_RUTI==0,31])),probs = c(0.25,0.75))
quantile(na.omit(as.numeric(table1_data[table1_RUTI==1,31])),probs = c(0.25,0.75))

## Table 1. 
## Patient characteristics related to UTI and RUTI (sample size = 963) used in the first stage analysis.
table1_result.m=rbind(age.m,gender.m,PC.m,dis.m,four_dis.m,
                      pre_hos.m,pre_uti_er.m,pre_uti_hos.m,
                      any_utisym.m,utisym.m,
                      creatinine.m,BWBC.m,UBact.m,UWBC.m,URBC.m)

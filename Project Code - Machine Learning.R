###PROJECT CODE###

####FINDING EPIGENETIC SIGNATURES FOR BIOLOGICAL AGEING####

#running the libraries required for the analysis of this project
library(MASS)
library(ggplot2)
library(GGally)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(kableExtra)
library(rpart)
library(rpart.plot)
library(caret)
library(randomForest)
library(glmnet)
library(forcats)
library(boot)
library(skimr)
library(ranger)
library(tidyverse)
library(caret)
library(pROC)
library(xgboost)

#Importing the dataset

data=read.csv("aDMPs_Proj2.csv",na.strings = NA)
head(data)
dim(data)[1]

#Checking for missing values
colSums(is.na(data))

#Prolif_meth and Senes_meth have missing values that need to be treated. Imputing these values using group means based on 
#aDMP status

#Treating the proliferating missing values
data = data%>%
  group_by(aDMP_status)%>%
  mutate(Prolif_Meth=ifelse(is.na(Prolif_Meth),mean(Prolif_Meth,na.rm=TRUE),Prolif_Meth))

#Checking with the prolif means
Prolif_means=data %>%
  group_by(aDMP_status)%>%
  summarise(Pmean=mean(Prolif_Meth,na.rm=TRUE))
Prolif_means
data[data$CpG_ID=="cg14008030",]

#Treating the senescent missing values
data = data%>%
  group_by(aDMP_status)%>%
  mutate(Senes_Meth=ifelse(is.na(Senes_Meth),mean(Senes_Meth,na.rm=TRUE),Senes_Meth))

#Checking with the senes means
Senes_means=data %>%
  group_by(aDMP_status)%>%
  summarise(Smean=mean(Senes_Meth,na.rm=TRUE))
Senes_means

#Rechecking if any missing values left
colSums(is.na(data))

#Changing the data into acceptable format for analysis
data$aDMP_status=as.factor(data$aDMP_status)
glimpse(data)

##Train/Test split using 70/30 ratio
set.seed(123)

n=nrow(data)
ind=sample(c(1:n),0.7*n)
data.train=data[ind,]
data.test=data[-ind,]

#Transforming the training dataset
constant=1
trans_data=data.train %>%
  mutate(across(c(2,15),~log(.x+constant))) #log transformation +1

#Scaling the training dataset
data.train.scale=trans_data%>%
  mutate(across(c(2,15),~(.-min(.))/max(.)))
glimpse(data.train.scale)

#Transforming the testing dataset
trans_data.test=data.test %>%
  mutate(across(c(2,15),~log(.x+constant))) #log transformation +1

#Scaling the testing dataset
data.test.scale=trans_data.test%>%
  mutate(across(c(2,15),~(.-min(.))/max(.)))
glimpse(data.test.scale)

###EXPLORATORY ANALYSIS

#Numerical summaries

#Statistics for proliferate cells 
#Summar for aDMP=0
my_skim <- skim_with(numeric = sfl(hist = NULL), 
                     base = sfl(n = length))

data.train.scale %>% 
  filter(aDMP_status==0)%>%
  my_skim() %>% 
  transmute(Variable=skim_variable, Mean=numeric.mean, SD=numeric.sd,
            Min=numeric.p0,  Max=numeric.p100,
            IQR = numeric.p75-numeric.p50) %>%
  kable(format.args = list(big.mark = ","),
        caption = '\\label{tab: Summary Statistics} 
        Summary statistics for Proliferate cells by aDMP Status', digits=4) %>%
  kable_styling(font_size = 10, latex_options = "hold_position")

#Summary for aDMP=1
data.train.scale %>% 
  filter(aDMP_status==1)%>%
  my_skim() %>% 
  transmute(Variable=skim_variable, Mean=numeric.mean, SD=numeric.sd,
            Min=numeric.p0,  Max=numeric.p100,
            IQR = numeric.p75-numeric.p50) %>%
  kable(format.args = list(big.mark = ","),
        caption = 'Summary statistics for aDMP sites', digits=4) %>%
  kable_styling(font_size = 10, latex_options = "hold_position")


#Graphical summaries
#Comparison based on aDMP_status - Density plots

Pplot_list=list() #Proliferate cells
Splot_list=list() #Senescent cells

P_data=data.train.scale[,2:8]
S_data=data.train.scale[,9:15]

for (i in 1:length(P_data)) {
  Pplot=ggplot(data.train.scale,aes(x=P_data[[i]], fill=aDMP_status))+
    geom_density()+
    labs(title=paste("Density Plot - ",names(P_data)[[i]]),x = names(P_data)[[i]])
  Pplot_list[[i]]=Pplot
  print(Pplot)
}

for (i in 1:length(S_data)) {
  Splot=ggplot(data.train.scale,aes(x=S_data[[i]],fill=aDMP_status))+
    geom_density()+
    labs(title=paste("Density Plot - ",names(S_data)[[i]]),x = names(S_data)[[i]])
  Splot_list[[i]]=Splot
  print(Splot)
}


#Comparison based on aDMP_status - Boxplots

PBplot_list=list() #Proliferate cells
SBplot_list=list() #Senescent cells

P_data=data.train.scale[,2:8]
S_data=data.train.scale[,9:15]

for (i in 1:length(P_data)) {
  PBplot=ggplot(data.train.scale,aes(y=P_data[[i]],x=aDMP_status))+
    geom_boxplot()+
    labs(title=paste("BoxPlot - Comparison of",names(P_data)[[i]], "based on aDMP status"),y = names(P_data)[[i]])
  PBplot_list[[i]]=PBplot
  print(PBplot)
}

for (i in 1:length(S_data)) {
  SBplot=ggplot(data.train.scale,aes(y=S_data[[i]], x=aDMP_status))+
    geom_boxplot()+
    labs(title=paste("BoxPlot - Comparison of",names(S_data)[[i]], "based on aDMP status"),y = names(S_data)[[i]])
  SBplot_list[[i]]=SBplot
  print(SBplot)
}

#Comparison based on Young or Old cells

PSplot_list=list()

for (i in 1:length(P_data)) {
  PSplot=ggplot(data.train.scale) +
    geom_density(aes(x=P_data[[i]],color="Proliferate"),alpha=0.5)+
    geom_density(aes(x=S_data[[i]],color="Senescent"),alpha=0.5)+
    labs(title = "Comparison of Proliferate and Senescent cells",x=paste(names(P_data)[[i]], "&" ,names(S_data)[[i]]), color="Cell Type")+
    scale_color_manual(values=c("red","blue"))
  PSplot_list[[i]]=PSplot
  print(PSplot)
}

#Pairsplots

Ppairs_plot=cbind(P_data,data.train.scale$aDMP_status)
Spairs_plot=cbind(S_data,data.train.scale$aDMP_status)

ggpairs(Ppairs_plot)
ggpairs(Spairs_plot)

#Lasso regression for feature selection
x=as.matrix(data.train.scale[,c(2:15)])
y=as.numeric(data.train.scale$aDMP_status)

lassomodel=cv.glmnet(x,y,family="gaussian",alpha=1)
optimal_lambda=lassomodel$lambda.min
selected_features=coef(lassomodel,s=optimal_lambda,exact=TRUE)
selected_features
feature_importance=rowSums(abs(selected_features))[-1]
ranked_features=sort(feature_importance,decreasing = TRUE)
ranked_features
plot(lassomodel)

#As per the lasso model results, all features have a non-zero coefficient. Thus, all features are significant as this point.
#We can do further feature selection based on our fitted classification models

###MODEL FITTING

###Logistic regression with bootstrapping

batch_size=30000
num_batches_lr=ceiling(nrow(data.train.scale)/batch_size)
log_models_lr <- list()

X=data.train.scale[,2:15]
Y=data.train.scale$aDMP_status

num_models_lr=ceiling(nrow(data.train.scale)/batch_size)
for (i in 1:num_models_lr) {
  bootstrapped_indices <- sample(1:nrow(data.train.scale), size = batch_size, replace = TRUE)
  bootstrapped_data <- data.train.scale[bootstrapped_indices,]
  X_batch <- bootstrapped_data[,2:15]
  Y_batch <- as.factor(bootstrapped_data$aDMP_status)
  log_model_lr <- glm(Y_batch~.,data=data.frame(X_batch,Y_batch),family=binomial)
  log_models_lr[[i]] <- log_model_lr
}
log_models_lr

# Step 3: Combining all logistic regression models
combine_models_lr <- function(models) {
  model_coeffs <- lapply(models, coef)
  avg_coeffs <- Reduce('+', model_coeffs) / length(models)
  final_model <- log_models_lr[[1]]
  final_model$coefficients <- avg_coeffs
  
  return(final_model)
}

final_model_lr <- combine_models_lr(log_models_lr)
final_model_lr
summary(final_model_lr)

#Odds ratio plot for all features

allcoefs=coef(final_model_lr)[-1]
allCI=confint(final_model_lr)[-1,]
allodds= data.frame(
  allvar = names(allcoefs),
  allodds_ratio = exp(allcoefs),
  lower_all = exp(allCI[, 1]),
  upper_all = exp(allCI[, 2]),
  alleffect=ifelse(allodds_ratio < 1, "Negative", "Positive")
)

ggplot(allodds, aes(x = allodds_ratio, y = allvar, color = alleffect)) +
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = lower_all, xmax = upper_all), width = 0.2) +
  labs(x = "Odds Ratio", y = "Variable",
       title = "Odds Ratio Plot") +
  scale_color_manual(values = c("Positive" = "green", "Negative" = "red")) +  # Color mapping
  theme_minimal()+
  scale_x_continuous(limits = c(0,range(allodds_ratio)[2]+0.1))

#The data seems to be highly imbalanced.
#Checking the ratio
ratio=sum(data.train.scale$aDMP_status==1)/(sum(data.train.scale$aDMP_status==0)+
                                              sum(data.train.scale$aDMP_status==1))
ratio*100
sum(data$aDMP_status==1)

#Predicting on train data
train_class_lr=as.factor(ifelse(train_predictions_lr>=optimal_threshold_lr,1,0))
head(train_class_lr)
confusionMatrix(train_class_lr,Y)
error_metrics(table(train_class_lr,Y))

#Predicting on test data
X_test=data.test.scale[,2:15]
Y_test=as.factor(data.test.scale$aDMP_status)
test_prediction_lr=predict(final_model_lr,X_test,type="response")
test_class_lr=as.factor(ifelse(test_prediction_lr>=optimal_threshold_lr,1,0))
confusionMatrix(test_class_lr,Y_test)

#ROC Curve
lr.rocCurve=roc(Y_test,test_prediction_lr)
plot(lr.rocCurve, main = "ROC Curve", col = "red", lwd = 2,xlab="False positive rate",ylab="True positive rate")
auc(lr.rocCurve)
#0.8684

###Logistic regression with class weights as data appears to be imbalanced

class_freq=table(data.train.scale$aDMP_status)
total_samp=sum(class_freq)
class_weights_cw <- total_samp / (length(class_freq) * class_freq)

batch_size=30000
num_batches_cw=ceiling(nrow(data.train.scale)/batch_size)
log_models_cw <- list()

num_models_cw=ceiling(nrow(data.train.scale)/batch_size)
for (i in 1:num_models_cw) {
  bootstrapped_indices <- sample(1:nrow(data.train.scale), size = batch_size, replace = TRUE)
  bootstrapped_data <- data.train.scale[bootstrapped_indices,]
  X_batch <- bootstrapped_data[,2:15]
  Y_batch <- as.factor(bootstrapped_data$aDMP_status)
  class_freq_cw <- table(Y_batch)
  total_samp_cw <- sum(class_freq_cw)
  class_weights_cw <- total_samp_cw / (length(class_freq_cw) * class_freq_cw)
  log_model_cw <- train(X_batch, Y_batch, method = "glm", family = "binomial", 
                        weights = ifelse(Y_batch == 1, class_weights_cw[2], class_weights_cw[1]))
  log_models_cw[[i]] <- log_model_cw
}
log_models_cw

# Step 3: Combining all logistic regression models
combine_models_cw <- function(models) {
  model_coeffs <- lapply(models, coef)
  avg_coeffs <- Reduce('+', model_coeffs) / length(models)
  final_model <- log_models_cw[[1]]
  final_model$coefficients <- avg_coeffs
  
  return(final_model)
}

final_model_cw <- combine_models_cw(log_models_cw)
summary(final_model_cw)

#ROC Curve
test_prediction_cw=predict(final_model_cw,X_test,type="prob")[,2]
lrcw.rocCurve=roc(Y_test,test_prediction_cw)
plot(lrcw.rocCurve, main = "ROC Curve", col = "red", lwd = 2,xlab ="False positive rate",ylab="True positive rate")
auc(lrcw.rocCurve)
#0.8613

###Logistic Regression with downsampling
X_down=data.train.scale[,2:15]
Y_down=data.train.scale$aDMP_status
downsample_data=downSample(X_down,Y_down)
head(downsample_data)
logdown_model=glm(Class~.,data=downsample_data,family=binomial)
summary(logdown_model)


#ROC Curve
test_prediction_ds=predict(logdown_model,downsample_data[,1:14],type="response")
lrds.rocCurve=roc(downsample_data$Class,test_prediction_ds)
plot(lrds.rocCurve, main = "ROC Curve", col = "red", lwd = 2,xlab="False positive rate",ylab="True positive rate")
auc(lrds.rocCurve)
#0.8985

#Since the AUC is a little low, but the AIC is too high, using the logistic regression model without class
#Weights or downsampling would be preferrable.

###Logistic Regression with feature selection using the original model

#Feature selection has to be done based on p-value. Selecting only those features with p-value less than 0.05.
#This gives us Prolif_H3.3, Prolif_Meth, Senes_H3.3,Senes_H4K16ac_ab1 and Senes_Meth.

final_data.train=data.train.scale[,c("Prolif_H3.3","Prolif_Meth","Senes_H3.3","Senes_H4K16ac_ab1","Senes_Meth","aDMP_status")]
final_data.test=data.test.scale[,c("Prolif_H3.3","Prolif_Meth","Senes_H3.3","Senes_H4K16ac_ab1","Senes_Meth","aDMP_status")]

batch_size=30000
num_batches_fs=ceiling(nrow(final_data.train)/batch_size)
log_models_fs <- list()

X_final=final_data.train[,1:5]
Y_final=final_data.train$aDMP_status

num_models_fs=ceiling(nrow(final_data.train)/batch_size)
for (i in 1:num_models_fs) {
  bootstrapped_indices_fs <- sample(1:nrow(final_data.train), size = batch_size, replace = TRUE)
  bootstrapped_data_fs <- final_data.train[bootstrapped_indices_fs,]
  X_batch_fs <- bootstrapped_data_fs[,1:5]
  Y_batch_fs <- as.factor(bootstrapped_data_fs$aDMP_status)
  log_model_fs <- glm(Y_batch_fs~.,data=data.frame(X_batch_fs,Y_batch_fs),family=binomial)
  log_models_fs[[i]] <- log_model_fs
}
log_models_fs

# Step 3: Combining all logistic regression models
combine_models_fs <- function(models) {
  model_coeffs <- lapply(models, coef)
  avg_coeffs <- Reduce('+', model_coeffs) / length(models)
  final_model <- log_models_fs[[1]]
  final_model$coefficients <- avg_coeffs
  
  return(final_model)
}

final_model_fs <- combine_models_fs(log_models_fs)
summary(final_model_fs)

#Plots for odds-ratios
coefs=coef(final_model_fs)[-1]
CI=confint(final_model_fs)[-1,]
odds= data.frame(
  variable = names(coefs),
  odds_ratio = exp(coefs),
  lower = exp(CI[, 1]),
  upper = exp(CI[, 2]),
  effect=ifelse(odds_ratio < 1, "Negative", "Positive")
)

ggplot(odds, aes(x = odds_ratio, y = variable, color = effect)) +
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2) +
  labs(x = "Odds Ratio", y = "Variable",
       title = "Odds Ratio Plot") +
  scale_color_manual(values = c("Positive" = "green", "Negative" = "red")) +  # Color mapping
  theme_minimal()+
  scale_x_continuous(limits = c(0,range(odds_ratio)[2]+0.1))

#ROC Curve
X_test_final=final_data.test[,1:5]
Y_test_final=as.factor(final_data.test$aDMP_status)
test_prediction_fs=predict(final_model_fs,X_test_final,type="response")
lrfs.rocCurve=roc(Y_test_final,test_prediction_fs)
plot(lrfs.rocCurve, main = "ROC Curve", col = "red", lwd = 2,xlab="False positive rate",ylab="True positive rate")
auc(lrfs.rocCurve)
#0.8687

###Random Forest using entire data, just to confirm if the selected features from logistic regression are 
#matching with the rf results

data.tree=data.train.scale[,c(2:16)]
data.tree.test=data.test.scale[,c(2:16)]

var_counts=c(3,4,5,6,7,8)
error=c()

for (i in 1:length(var_counts)) {
  bag_model=randomForest(aDMP_status~.,data = data.tree,mtry=var_counts[i],ntree=300)
  error[i]=mean(predict(bag_model,data.tree,type="class")!=data.tree$aDMP_status)
}
error
var.opt=which.min(error)
var.opt

bag_model.opt= randomForest(aDMP_status~.,data = data.tree,mtry=var.opt,ntree=300)

bag.train=predict(bag_model.opt,data.tree,type="class")
bag.test=predict(bag_model,data.tree.test,type="class")
#ROC Curve
rf.probs=predict(bag_model.opt,data.tree.test,type="prob")
rf.rocCurve=roc(data.tree.test$aDMP_status,rf.probs[,2])
plot(rf.rocCurve,col=c(2),main = "ROC Curve",xlab="False positive rate",ylab="True positive rate")
auc(rf.rocCurve)
#0.8623

#Feature selection
bag_df=data_frame(var=rownames(randomForest::importance(bag_model.opt)),MeanDecreaseGini=randomForest::importance(bag_model.opt)[,1])%>%
  mutate(var=fct_reorder(var,MeanDecreaseGini,median))
bag_ggplot=ggplot(bag_df,aes(var,MeanDecreaseGini))+
  geom_point()+
  coord_flip()+
  labs(title="Feature selection through Gini Index")
bag_ggplot
sort(bag_df$MeanDecreaseGini,decreasing = TRUE)

#The features selected from logistic regression model are similar to the top features from rf model

###Random forest with selected features
var_counts_rf=c(2,3,4,5)
error_rf=c()

for (i in 1:length(var_counts_rf)) {
  bag_model_rf=randomForest(aDMP_status~.,data = final_data.train,mtry=var_counts_rf[i],ntree=300)
  error_rf[i]=mean(predict(bag_model_rf,final_data.train,type="class")!=final_data.train$aDMP_status)
}

var.opt.rf=which.min(error_rf)
var.opt.rf

bag_model.final= randomForest(aDMP_status~.,data = final_data.train,mtry=var.opt.rf,ntree=300) 

#ROC Curve
set.seed(123)
rf.probs.train=predict(bag_model.final,final_data.train,type="prob")
rf.train_rocCurve=roc(final_data.train$aDMP_status,rf.probs.train[,2])
plot(rf.train_rocCurve,col=c(2),main="ROC Curve",xlab="False positive rate",ylab="True positive rate")
auc(rf.train_rocCurve)
#0.9986

rf.probs.test=predict(bag_model.final,final_data.test,type="prob")
rf.test_rocCurve=roc(final_data.test$aDMP_status,rf.probs.test[,2])
plot(rf.test_rocCurve,col=c(2),main="ROC Curve",xlab="False positive rate",ylab="True positive rate")
auc(rf.test_rocCurve)
#0.8432

###Random forest with class weights
class_freq=table(final_data.train$aDMP_status)
total_samp=sum(class_freq)
class_weights <- total_samp / (length(class_freq) * class_freq)

error_rfcw=c()
for (i in 1:length(var_counts_rf)) {
  bag_model_cw=randomForest(aDMP_status~.,data =final_data.train,mtry=var_counts_rf[i],ntree=300,
                            classwt=class_weights)
  error_rfcw[i]=mean(predict(bag_model_cw,final_data.train,type="class")!=final_data.train$aDMP_status)
}

var.opt.rfcw=which.min(error_rfcw)
var.opt.rfcw

bag_model_cw.opt= randomForest(aDMP_status~.,data = final_data.train,mtry=var.opt.rfcw,ntree=300,
                               classwt=class_weights)

#ROC Curve
set.seed(123)
rfcw.probs=predict(bag_model_cw.opt,final_data.train,type="prob")
rfcw_rocCurve=roc(final_data.train$aDMP_status,rfcw.probs[,2])
plot(rfcw_rocCurve,col=c(2),main="ROC Curve",xlab="False positive rate",ylab="True positive rate")
auc(rfcw_rocCurve)
#0.9996

rfcw.probs.test=predict(bag_model_cw.opt,final_data.test,type="prob")
rfcw.test_rocCurve=roc(final_data.test$aDMP_status,rfcw.probs.test[,2])
plot(rfcw.test_rocCurve,col=c(2),main="ROC Curve",xlab="False positive rate",ylab="True positive rate")
auc(rfcw.test_rocCurve)
#0.842


###Gradient Boosting
set.seed(123)

model.gbm=train(aDMP_status~.,data=final_data.train,method="xgbTree", trControl=trainControl("cv",number=10))

#ROC Curve

gbm.probs=predict(model.gbm,final_data.test,type="prob")
head(rf.probs)

gbm.rocCurve=roc(final_data.test$aDMP_status,gbm.probs[,2])
plot(gbm.rocCurve,col=c(2),main="Area under the curve: 0.911",xlab="False positive rate",ylab="True positive rate")
auc(gbm.rocCurve)
#0.9108

#Feature importance
trained_model <- model.gbm$finalModel
importance_values <- xgb.importance(model = trained_model)

xgb.plot.importance(importance_matrix = importance_values,main="Feature importance")

###Gradient Boosting with class weights

set.seed(123)
class_freq_gbm=table(final_data.train$aDMP_status)
total_samp_gbm=sum(class_freq_gbm)
class_weights_gbm <- total_samp_gbm / (length(class_freq_gbm) * class_freq_gbm)

model.gbm_cw=train(aDMP_status~.,data=final_data.train,method="xgbTree", 
                   trControl=trainControl("cv",number=10),
                   weights = ifelse(final_data.train$aDMP_status == 1, 
                                    class_weights_gbm[2], class_weights_gbm[1]))

#ROC Curve

gbm_cw.probs=predict(model.gbm_cw,final_data.test,type="prob")

gbm_cw.rocCurve=roc(final_data.test$aDMP_status,gbm_cw.probs[,2])
plot(gbm_cw.rocCurve,col=c(2),main="ROC Curve",xlab="False positive rate",ylab="True positive rate")
auc(gbm_cw.rocCurve)
#0.9078

###Gradient boosting with hyperparameter tuning
param_grid <- expand.grid(
  n.trees = c(100, 200, 300),
  interaction.depth = c(2, 4, 8),
  shrinkage = c(0.01, 0.1, 0.4),
  n.minobsinnode = c(10000, 20000, 30000)
)

ctrl <- trainControl(
  method = "cv",          
  number = 10,            
  search = "grid",        
  verboseIter = TRUE
)

model_gbm_cv <- train(
  aDMP_status ~ .,                  
  data = final_data.train,       
  method = "gbm",         
  trControl = ctrl,       
  tuneGrid = param_grid   
)

gbm_cv.probs=predict(model_gbm_cv,final_data.test,type="prob")

gbm_cv.rocCurve=roc(final_data.test$aDMP_status,gbm_cv.probs[,2])
plot(gbm_cv.rocCurve,col=c(2),xlab = "False positive rate",ylab="True positive rate")
auc(gbm_cv.rocCurve)
#0.8891

###The best model out of these was the gradient boosting model with 10-fold cross validation.

#########Secondary aim

#Proliferating cells

Prolif_data.train=data.train.scale[,c(2:8,16)]
Prolif_data.test=data.test.scale[,c(2:8,16)]

#Proliferating cells using Gradient boosting

#Entire proliferate data
set.seed(123)

model.gbm.pfull=train(aDMP_status~.,data=Prolif_data.train,method="xgbTree", trControl=trainControl("cv",number=10))

#Feature importance
Prolif_fullmodel <- model.gbm.pfull$finalModel
feat_imp_pfull <- xgb.importance(model = Prolif_fullmodel)

xgb.plot.importance(importance_matrix = feat_imp_pfull,main="Feature importance for proliferate cells")

#Selected features
PFinal_data.train=Prolif_data.train[,c("Prolif_H3.3","Prolif_Meth","aDMP_status")]
PFinal_data.test=Prolif_data.test[,c("Prolif_H3.3","Prolif_Meth","aDMP_status")]

set.seed(123)

model.gbm.pf=train(aDMP_status~.,data=PFinal_data.train,method="xgbTree", trControl=trainControl("cv",number=10))

#ROC Curve
gbm.pf.probs=predict(model.gbm.pf,Prolif_data.test,type="prob")

gbm.pf.rocCurve=roc(Prolif_data.test$aDMP_status,gbm.pf.probs[,2])
plot(gbm.pf.rocCurve,col=c(2),main="AUC is 0.846 for proliferate cells",ylab="True positive rate",
     xlab="True negative rate")
auc(gbm.pf.rocCurve)
#0.8456

#Senescent cells

Senes_data.train=data.train.scale[,9:16]
Senes_data.test=data.test.scale[,9:16]

#Senescent cells using Gradient boosting
set.seed(123)
#entire data
model.gbm.sfull=train(aDMP_status~.,data=Senes_data.train,method="xgbTree", trControl=trainControl("cv",number=10))

#Feature importance
Senes_fullmodel <- model.gbm.sfull$finalModel
feat_imp_sfull <- xgb.importance(model = Senes_fullmodel)

xgb.plot.importance(importance_matrix = feat_imp_sfull,main="Feature importance for senescent cells")

#Selected features

SFinal_data.train=Senes_data.train[,c("Senes_H3.3","Senes_H4K20me3_ab1","Senes_Meth","aDMP_status")]
SFinal_data.test=Senes_data.test[,c("Senes_H3.3","Senes_H4K20me3_ab1","Senes_Meth","aDMP_status")]

set.seed(123)

model.gbm.sn=train(aDMP_status~.,data=SFinal_data.train,method="xgbTree", trControl=trainControl("cv",number=10))

#ROC Curve
gbm.sn.probs=predict(model.gbm.sn,Senes_data.test,type="prob")

gbm.sn.rocCurve=roc(Senes_data.test$aDMP_status,gbm.sn.probs[,2])
plot(gbm.sn.rocCurve,col=c(2),main="AUC is 0.679 for senescent cells",ylab="True positive rate",
     xlab="True negative rate")
auc(gbm.sn.rocCurve)
#0.6873


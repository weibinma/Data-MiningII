#install.packages("psych")
rm(list = ls())
library(ggplot2)
library(psych)
library(ElemStatLearn)
library(MASS)
library(arules)
source("http://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("multtest")
#biocLite("cluster")
#install.packages("fpc")
library("multtest")
library("fpc")
library("cluster")
#library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(multtest)
# load the original data
data <- read.csv("output.csv")
#head(data)
dataset <- data[,5:84]
#check scaling
layout(matrix(c(1,1),nrow=1, ncol=2))
boxplot(dataset, main= "boxplot for all 80 samples")
dim(dataset)
#transpose dataset
dataset_t = t(dataset)
dim(dataset_t)
#doing pca
datapca <- prcomp(dataset_t,center=TRUE, t=0.01)
dim(datapca$x)
dim(datapca$rotation)
summary(datapca$x)
quartz()
plot(datapca)

#choose the optimal value of k
##k-means_gap
gap_km <- clusGap(datapca$x, kmeans, nstart = 1, K.max = 10, B = 100)
quartz()
plot(gap_km, main = "Gap Statistic: k-means")
##k-medoids_gap
gap_kmed <- clusGap(datapca$x, pam, K.max = 10, B = 100)
quartz()
plot(gap_kmed, main = "Gap Statistic: k-medoids")
#k-means_wss
set.seed(123)
quartz()
fviz_nbclust(datapca$x, kmeans, method = "wss")
# x11()
# fviz_nbclust(datapca$x, kmeans, method = "silhouette")

# Kmeans_k=4
distance <- dist(datapca$x)
km <- kmeans(datapca$x, 4, nstart = 20)
si_kmeans <- silhouette(km$cluster, dist = distance)
newtype=cbind(dataset_t[,0],km$cluster)
summary(si_kmeans)
quartz()
plot(si_kmeans,col=2:5)
tk<- table(km$cluster)
tk
quartz()
fviz_cluster(km, data = datapca$x)

# Kmeans_k=5
km1 <- kmeans(datapca$x, 5, nstart = 20)
si_kmeans1 <- silhouette(km1$cluster, dist = distance)
k_5 <- dataset_t[,0]
newtype1 = cbind(dataset_t[,0], km1$cluster)
newtype1 = as.data.frame(km1$cluster)

summary(si_kmeans1)
quartz()
plot(si_kmeans1,col=2:6)
tk1<- table(km1$cluster)
tk1
quartz()
fviz_cluster(km1, data = datapca$x)




##########################################################################################################



dataset_proteomes <- read.csv('77_cancer_proteomes_CPTAC_itraq.csv')
dataset_clinical <- read.csv('clinical_data_breast_cancer.csv')
#...decide which Tumor is the most fatal leading to Breast Cancer. Here I suppose Tumor2 as the most fatal one.
#dataset_T2 <- dataset[which(dataset$Tumor=='T2'),]

###########################
#find the most fatal tumor.
###########################
#remove unimportant columns
dataset_survived <- dataset_clinical[,c(4,8,10,11,15,17,22)]
#rename 'Age', 'AJCC Status' column
colnames(dataset_survived)[1] <- 'Age'
colnames(dataset_survived)[6] <- 'Survived'

#set class for 'Age' column
dataset_survived$Age <- ordered(cut(dataset_survived$Age, c(29,40,60,90), labels = c("young", "old", "older")))
#convert data into transactions
dataset_survived_tran <- as(dataset_survived, "transactions")
#summary(dataset_survived_tran)
#find item frequency
itemFrequency(dataset_survived_tran)
#visualize data
quartz()
itemFrequencyPlot(dataset_survived_tran, support = 0.1, ylim = range(0,1))
#apply the apriori algorithm
apriori(dataset_survived_tran)
survived_tran_rules <- apriori(dataset_survived_tran, parameter = list(support = 0.016,
                                                             confidence = 0.85))
#find area schooling a priority and with low pupil-teacher ratios
dead_rules <- subset(survived_tran_rules, subset = rhs %in% "Survived=DECEASED" &lift>1)
dead_rules 
inspect(head(dead_rules))
##Conclusion: 
#1. patients' Node.Coded = Positive is more likely to death.
#2. older patients with Stage IV are more likely to death.
#3. Category 'Basal-like' is more likely to be the most fatal tumor in these four tumor categories


###########################
#the most fatal tumor.
###########################
#extract PAM50.mRNA=='Basal-like'
targeted_tumor <- dataset_clinical[which(dataset_clinical$PAM50.mRNA=='Basal-like'),]
#extract names of patients with 'Basal-like' tumor
#targeted_tumor_names <- substring(targeted_tumor$New.tumor.clusters, 1, 7)
#targeted_tumor_names <- paste('TCGA-',targeted_tumor_names, sep = '')

#####################
##supervised learning
#####################

#remove unimportant columns
dataset_survived <- dataset_clinical[,c(4,8,10,11,15,17,22)]
#rename 'Age', 'AJCC Status' column
colnames(dataset_survived)[1] <- 'Age'
colnames(dataset_survived)[4] <- 'Spread_Status'
colnames(dataset_survived)[6] <- 'Survived'

####predicting the breast cancer recurrence chance
dataset_predict <- dataset_survived
#reordering label to the last column
dataset_predict <- dataset_predict[, c(1,2,3,5,6,7,4)]
# Encoding categorical variables
dataset_predict$Spread_Status = factor(dataset_predict$Spread_Status,
                         levels = c('Positive', 'Negative'),
                         labels = c(1, 0))

dataset_predict$Tumor = factor(dataset_predict$Tumor,
                                  levels = c('T1', 'T2', 'T3', 'T4'),
                                  labels = c(1,2,3,4))

dataset_predict$Node = factor(dataset_predict$Node,
                                  levels = c('N0', 'N1', 'N2', 'N3'),
                                  labels = c(1,2,3,4))

dataset_predict$Converted.Stage = factor(dataset_predict$Converted.Stage,
                                  levels = c('Stage I', 'Stage IA', 'Stage IB', 'Stage II', 'Stage IIA', 'Stage IIB',
                                             'Stage III', 'Stage IIIA', 'Stage IIIB', 'Stage IIIC', 'Stage IV'),
                                  labels = c(1,2,3,4,5,6,7,8,9,10,11))

dataset_predict$Survived = factor(dataset_predict$Survived,
                                  levels = c('DECEASED', 'LIVING'),
                                  labels = c(0,1))

dataset_predict$PAM50.mRNA = factor(dataset_predict$PAM50.mRNA,
                                  levels = c('Basal-like', 'HER2-enriched', 'Luminal A', 'Luminal B'),
                                  labels = c(1,2,3,4))

#scale
dataset_predict[,1:6] <- scale(data.matrix(dataset_predict[,1:6]))  #country and purchased entries is defined in "factor", not numeric

###### KNN
library(caTools)
set.seed(123)
split = sample.split(dataset_predict$Recurrence_label, SplitRatio = 0.8)
training_set = subset(dataset_predict, split == TRUE)
test_set = subset(dataset_predict, split == FALSE)

# Fitting K-NN to the Training set and Predicting the Test set results
library(class)
y_pred_knn = knn(train = training_set[, -7],
             test = test_set[, -7],
             cl = training_set[, 7],
             k = 5,
             prob = TRUE)

# Making the Confusion Matrix
cm_knn = table(test_set[, 7], y_pred_knn)
cm_knn #100%
confusionMatrix(cm_knn)

###### Decision Tree
set.seed(456)
library(randomForest)
classifier = randomForest(x = training_set[-7],
                          y = training_set$Recurrence_label,
                          ntree = 500)

# Predicting the Test set results
y_pred_decisionTree = predict(classifier, newdata = test_set[-7])

# Making the Confusion Matrix
cm_deci = table(test_set[, 7], y_pred_decisionTree)
cm_deci #100%
confusionMatrix(cm_deci)

###### Logistic Regression
# Fitting Logistic Regression to the Training set
classifier = glm(formula = Recurrence_label ~ .,
                 family = binomial,
                 data = training_set)

# Predicting the Test set results
prob_pred = predict(classifier, type = 'response', newdata = test_set[-7])
y_pred_log = ifelse(prob_pred > 0.5, 1, 0)

# Making the Confusion Matrix
cm_log = table(test_set[, 7], y_pred_log > 0.5)
cm_log #100%



#########################################################################################################

#rename
k_5 <- dataset_t[,0:1]
newtype_k_5 <- substring(names(k_5), 1, 7)
newtype_k_5 <- paste('TCGA-',newtype_k_5, sep = '')
substring(newtype_k_5, 8) <- '-'
newtype_k_5

#transpose gene data 
dataset_proteomes_t <- as.data.frame(dataset_t)
#save(dataset_proteomes_t, file="dataset_proteomes_t.csv")
dataset_proteomes_t['Complete.TCGA.ID'] <- newtype_k_5

#merge
df1 <- dataset_proteomes_t
df2 <- dataset_clinical[,c(2, 22)]
merge_dataset <- merge(df1, df2, by = 'Complete.TCGA.ID')
write.csv(merge_dataset, 'merge_80and105.csv')
############################
##modeling classification
############################
dim(merge_dataset)
######SVM without using dimension reduction
# Load library
library(dplyr)
library(caTools)
library(e1071)
library(caret) 
# Transform dataset
dataset_svm <- merge_dataset[, -1]
dataset_svm$PAM50.mRNA = factor(dataset_svm$PAM50.mRNA,
                               levels = c('Basal-like', 'Luminal A', 'Luminal B', 'HER2-enriched'),
                               labels = c(1,2,3,4))
# Splitting the dataset into the Training set and Test set
set.seed(789)
split = sample.split(dataset_svm$PAM50.mRNA, SplitRatio = 0.75)
training_set = subset(dataset_svm, split == TRUE)
test_set = subset(dataset_svm, split == FALSE)
#train data
classifier = svm(formula = PAM50.mRNA ~ .,
                 data = training_set,
                 type = 'C-classification',
                 kernel = 'linear')
# Predicting the Test set results
y_pred = predict(classifier, newdata = test_set[-9712])
# Making the Confusion Matrix
cm_svm = table(test_set[, 9712], y_pred)
cm_svm
confusionMatrix(cm_svm)  #65% accuracy
##Since only 65% accuracy without using dimension reduction, we decide to use LDA
##to reduce data dimension.

######SVM with using dimension reduction (LDA)
library(caTools)
set.seed(222)
#LDA
dataset_lda <- dataset_svm
split = sample.split(dataset_lda$PAM50.mRNA, SplitRatio = 0.8)
training_set_lda = subset(dataset_lda, split == TRUE)
test_set_lda = subset(dataset_lda, split == FALSE)

# Feature Scaling
training_set_lda[-9712] = scale(training_set_lda[-9712])
test_set_lda[-9712] = scale(test_set_lda[-9712])

# Applying LDA
library(MASS)
lda = lda(formula = PAM50.mRNA ~ ., data = training_set_lda)
training_set_lda = as.data.frame(predict(lda, training_set_lda))
training_set_lda = training_set_lda[c(6, 7, 8, 1)]
test_set_lda = as.data.frame(predict(lda, test_set_lda))
test_set_lda = test_set_lda[c(6, 7, 8, 1)]

# Fitting SVM to the Training set
# install.packages('e1071')
library(e1071)
classifier = svm(formula = class ~ .,
                 data = training_set_lda,
                 type = 'C-classification',
                 kernel = 'linear')

# Predicting the Test set results
y_pred_lda = predict(classifier, newdata = test_set_lda[-4])

# Making the Confusion Matrix
cm_svm_lda = table(test_set_lda[, 4], y_pred_lda)
cm_svm_lda
confusionMatrix(cm_svm_lda) #88.2% accuracy


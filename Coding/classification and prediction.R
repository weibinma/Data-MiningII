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

######loading dataset
# load the original data
data <- read.csv("output.csv")
dataset_clinical <- read.csv('clinical_data_breast_cancer.csv')
#head(data)
dataset <- data[,5:84]
#check scaling
layout(matrix(c(1,1),nrow=1, ncol=2))
boxplot(dataset, main= "boxplot for all 80 samples")
dim(dataset)
#transpose dataset
dataset_t = t(dataset)
dim(dataset_t)

#########################################################################################################

#rename
k_5 <- dataset_t[,0:1]
newtype_k_5 <- substring(names(k_5), 1, 7)
newtype_k_5 <- paste('TCGA-',newtype_k_5, sep = '')
substring(newtype_k_5, 8) <- '-'
newtype_k_5

#obtaining gene data
dataset_77cancer <- read.csv('output.csv')
dataset_77cancer <- dataset_77cancer[,-1, drop=FALSE]  #remove 1st column
dataset_77cancer <- dataset_77cancer[,-(2:3), drop=FALSE]  #remove 2,3 column
dataset_77cancer <- dataset_77cancer[,-(82:84), drop=FALSE]  #remove 82,83,84 column
colnames(dataset_77cancer)[1] <- 'GeneID'
colnames(dataset_77cancer)
dim(dataset_77cancer) #9711 81 ,which including 'GeneID'

#rename patients
k_original <- names(dataset_77cancer)[2:81]
k_original <- substring(k_original, 1, 7)
k_original <- paste('TCGA-',k_original, sep = '')
substring(k_original, 8) <- '-'
#k_original
names(dataset_77cancer)[2:81] <- k_original
names(dataset_77cancer) #9711 81 ,which including 'GeneID'
#transpose 
#dataset_77cancer_t <- t(dataset_77cancer)
#dataset_proteomes_t <- as.data.frame(dataset_77cancer_t)
#rownames(dataset_proteomes_t)
#dim(dataset_proteomes_t) #81 9711, which including 'GeneID'

#obtaining pam50 data
dataset_pam50 <- read.csv('PAM50_proteins.csv')
colnames(dataset_pam50)[2] <- 'GeneID'
names(dataset_pam50) #"GeneSymbol" "GeneID"     "Species"    "Gene.Name" 
dim(dataset_pam50) #100 4

#merge
df1 <- dataset_77cancer
df2 <- dataset_pam50[2]  #如果要加class，在这里加一列class就可以
original_merge_dataset <- merge(df1, df2, by = 'GeneID')      # 35 81 which means only 35 common gene.
write.csv(original_merge_dataset, 'original_merge_dataset.csv')

############################
##modeling classification
############################
dim(original_merge_dataset) #35 81
original_merge_dataset <- read.csv('original_merge_dataset.csv')
#transpose
original_merge_dataset_t <- as.data.frame(t(original_merge_dataset))
original_merge_dataset_t <- original_merge_dataset_t[-1,, drop=FALSE]  #remove 1st column
######SVM without using dimension reduction
# Load library
library(dplyr)
library(caTools)
library(e1071)
library(caret) 
# Transform dataset
dataset_svm <- read.csv('merge_dataset.csv')   #80 9714
dataset_svm <- dataset_svm[, -c(1,2)] #80 9712, which 9711 genes + 1 label
dataset_svm$PAM50.mRNA = factor(dataset_svm$PAM50.mRNA,
                                levels = c('Basal-like', 'Luminal A', 'Luminal B', 'HER2-enriched'),
                                labels = c(1,2,3,4))

######SVM without using dimension reduction (LDA)
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
a_svm <- confusionMatrix(cm_svm)
svm_acc <- a_svm$overall[1]  #65% accuracy

######Random Forest without using dimension reduction (LDA)
set.seed(123)
split = sample.split(dataset_svm$PAM50.mRNA, SplitRatio = 0.75)
training_set = subset(dataset_svm, split == TRUE)
test_set = subset(dataset_svm, split == FALSE)
# Feature Scaling
training_set[-9712] = scale(training_set[-9712])
test_set[-9712] = scale(test_set[-9712])
# Fitting Random Forest Classification to the Training set
# install.packages('randomForest')
library(randomForest)
set.seed(123)
classifier = randomForest(x = training_set[-9712],
                          y = training_set$PAM50.mRNA,
                          ntree = 500)
# Predicting the Test set results
y_pred_randf = predict(classifier, newdata = test_set[-9712])
# Making the Confusion Matrix
cm_randf = table(test_set[, 9712], y_pred_randf)
a_randf <- confusionMatrix(cm_randf)
randomForest_acc <- a_randf$overall[1] #0.7 accuracy

######KNN without using dimension reduction (LDA)
# Fitting K-NN to the Training set and Predicting the Test set results
library(class)
y_pred_knn = knn(train = training_set[, -9712],
             test = test_set[, -9712],
             cl = training_set[, 9712],
             k = 5,
             prob = TRUE)

# Making the Confusion Matrix
cm_knn = table(test_set[, 9712], y_pred_knn)
a_knn <- confusionMatrix(cm_knn)
knn_acc <- a_knn$overall[1] #0.45 accuracy

######Naive Bayes without using dimension reduction (LDA)
library(e1071)
classifier = naiveBayes(x = training_set[-9712],
                        y = training_set$PAM50.mRNA)
# Predicting the Test set results
y_pred = predict(classifier, newdata = test_set[-9712])
# Making the Confusion Matrix
cm_nai = table(test_set[, 9712], y_pred)
a_nai <- confusionMatrix(cm_nai)
naiveBayes_acc <- a_nai$overall[1]  #0.55 accuracy
 

######SVM with using dimension reduction (LDA)
library(caTools)
set.seed(789)
#LDA
dataset_lda <- dataset_svm
split = sample.split(dataset_lda$PAM50.mRNA, SplitRatio = 0.75)
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



######SVM with using dimension reduction (PCA)
library(caTools)
set.seed(789)
#PCA
dataset_pca <- dataset_svm
split = sample.split(dataset_pca$PAM50.mRNA, SplitRatio = 0.75)
training_set_pca = subset(dataset_pca, split == TRUE)
test_set_pca = subset(dataset_pca, split == FALSE)
# Feature Scaling
training_set_pca[-9712] = scale(training_set_pca[-9712])
test_set_pca[-9712] = scale(test_set_pca[-9712])
# Applying PCA
library(caret)
library(e1071)
pca = preProcess(x = training_set_pca[-9712], method = 'pca', pcaComp = 60)
training_set_pca = predict(pca, training_set_pca)
test_set_pca = predict(pca, test_set_pca)



########## Fitting SVM with lda
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
a_svm_lda <- confusionMatrix(cm_svm_lda)
svm_lda_acc <- a_svm_lda$overall[1] #85% accuracy

########## Fitting SVM with pca
classifier_svmpca = svm(formula = PAM50.mRNA ~ .,
                 data = training_set_pca,
                 type = 'C-classification',
                 kernel = 'linear')

# Predicting the Test set results
y_pred_svmpca = predict(classifier_svmpca, newdata = test_set_pca[-1])

# Making the Confusion Matrix
cm_svmpca = table(test_set_pca[, 1], y_pred_svmpca)
svmpca <- confusionMatrix(cm_svmpca)
svmpca #0.15


########## Fitting random forest with lda
set.seed(123)
classifier = randomForest(x = training_set_lda[,-4],
                          y = training_set_lda[,4],
                          ntree = 500)
# Predicting the Test set results
y_pred_randf = predict(classifier, newdata = test_set_lda[-4])
# Making the Confusion Matrix
cm_randf = table(test_set_lda[,4], y_pred_randf)
a_ranf_lda <- confusionMatrix(cm_randf)
randomForest_lda_acc <- a_ranf_lda$overall[1] #0.9 accuracy

########## Fitting random forest with pca
classifier_rfpca = randomForest(x = training_set_pca[,-1],
                          y = training_set_pca[,1],
                          ntree = 500)
# Predicting the Test set results
y_pred_rfpca = predict(classifier_rfpca, newdata = test_set_pca[-1])
# Making the Confusion Matrix
cm_rfpca = table(test_set_pca[,1], y_pred_rfpca)
rfpca <- confusionMatrix(cm_rfpca)
rfpca  #0.45


########## Fitting KNN with lda
y_pred_knn1 = knn(train = training_set_lda[, -4],
                 test = test_set_lda[, -4],
                 cl = training_set_lda[, 4],
                 k = 5,
                 prob = TRUE)
# Making the Confusion Matrix
cm_knn_new = table(test_set_lda[, 4], y_pred_knn1)
a_knn_lda <- confusionMatrix(cm_knn_new)
knn_lda_acc <- a_knn_lda$overall[1] #0.85 accuracy

########## Fitting KNN with pca
y_pred_knn2 = knn(train = training_set_pca[, -1],
                  test = test_set_pca[, -1],
                  cl = training_set_pca[, 1],
                  k = 5,
                  prob = TRUE)
# Making the Confusion Matrix
cm_knnpca = table(test_set_pca[, 1], y_pred_knn2)
knnpca <- confusionMatrix(cm_knnpca)
knnpca #0.75



########## Fitting Naive Bayes with lda
classifier = naiveBayes(x = training_set_lda[, -4],
                        y = training_set_lda[, 4])
# Predicting the Test set results
y_pred_nai = predict(classifier, newdata = test_set_lda[-4])
# Making the Confusion Matrix
cm_nai1 = table(test_set_lda[, 4], y_pred_nai)
a_nai_lda <- confusionMatrix(cm_nai1)
naiveBayes_lda_acc <- a_nai_lda$overall[1] #0.85
########## Fitting Naive Bayes with pca
classifierpca = naiveBayes(x = training_set_pca[, -1],
                        y = training_set_pca[, 1])
# Predicting the Test set results
y_pred_naipca = predict(classifierpca, newdata = test_set_pca[-1])
# Making the Confusion Matrix
naive_pca = table(test_set_pca[, 1], y_pred_naipca)
naipca <- confusionMatrix(naive_pca)
naipca #0.5



ro1 <- cbind(svm_acc, randomForest_acc, knn_acc, naiveBayes_acc)
ro2 <- cbind(svm_lda_acc, randomForest_lda_acc, knn_lda_acc, naiveBayes_lda_acc)
ro3 <- cbind(svmpca$overall[1], rfpca$overall[1], knnpca$overall[1], naipca$overall[1])
table_all <- rbind(ro1, ro2, ro3)
colnames(table_all) <- c('SVM', 'Random Forest', 'KNN', 'Naive Bayes')
rownames(table_all) <- c('Accuracy (No Dimension Reduction)', 'Accuracy (LDA)',
                         'Accuracy (PCA)')
table_all





######Visualizing
#install.packages("scatterplot3d") # Install
library("scatterplot3d") # load
colors <- c("#999999", "#E69F00", "#56B4E9", 'red')
all_points <- rbind.data.frame(training_set_lda, test_set_lda)
colors <- colors[as.numeric(all_points$class)]
scatterplot3d(all_points[,1:3], pch = 16, color=colors, grid=TRUE, box=FALSE)


#pca plot
all_points_pca <- rbind.data.frame(training_set_pca, test_set_pca)
colors <- c("#999999", "#E69F00", "#56B4E9", 'red')
colors <- colors[as.numeric(all_points_pca$PAM50.mRNA)]
scatterplot3d(all_points_pca[2:4], pch = 16, color=colors, grid=TRUE, box=FALSE)



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
names(dataset_77cancer)
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
original_merge_dataset <- merge(df1, df2, by = 'GeneID')
write.csv(original_merge_dataset, 'original_merge_dataset.csv')

############################
##modeling classification
############################
dim(original_merge_dataset) #35 81
#transpose
original_merge_dataset_t <- as.data.frame(t(original_merge_dataset))
original_merge_dataset_t <- original_merge_dataset_t[-1,, drop=FALSE]  #remove 1st column
colnames(original_merge_dataset_t) <- original_merge_dataset$GeneID
######SVM without using dimension reduction
# Load library
library(dplyr)
library(caTools)
library(e1071)
library(caret) 

######clustering
clu_dataset <- original_merge_dataset_t  #80 35
clu_dataset <- add_rownames(clu_dataset, 'Complete TCGA ID')
clu_dataset <- as.data.frame(clu_dataset)
new_clu_dataset <-clu_dataset[,-1]
new_clu_dataset = as.data.frame(sapply(new_clu_dataset, as.numeric))
new_clu_dataset <- scale(new_clu_dataset)
new_clu_dataset <- as_data_frame(new_clu_dataset)
#new_clu_dataset <- new_clu_dataset[-c(78,79,80),]
#choose the optimal value of k
set.seed(0428)
#wss
fviz_nbclust(new_clu_dataset, kmeans, method = "wss")  # k = 3,4,5
#sihouette
fviz_nbclust(new_clu_dataset, kmeans, method = "silhouette")  # k = 2,3,4
#gap
gap_km <- clusGap(new_clu_dataset, kmeans, nstart = 20, K.max = 10, B = 100)
print(gap_km, method = 'firstmax')
plot(gap_km, main = "Gap Statistic: k-means")
fviz_gap_stat(gap_km) # k = 3,4,5
#indices
library('NbClust')
##Description
#NbClust package provides 30 indices for determining the number of clusters and proposes
#to user the best clustering scheme from the different results obtained by varying all 
#combinations of number of clusters, distance measures, and clustering methods.
nb <- NbClust(new_clu_dataset, distance = "euclidean", min.nc = 2,
              max.nc = 10, method = "kmeans", index ="all")
fviz_nbclust(nb) + theme_minimal()  #10 proposed k = 3 is optimal, 7 proposed k = 2 is optimal.

#Thus, k = 3,4,5 might be optimal.

#####
#k=3
#####
distance <- dist(new_clu_dataset)
km_3 <- kmeans(new_clu_dataset, 3, nstart = 20)
si_kmeans <- silhouette(km_3$cluster, dist = distance)
summary(si_kmeans)
#silhouette plot
plot(si_kmeans, col=2:4)
tk_3 <- table(km_3$cluster)
#cluster plot
fviz_cluster(km_3, data = new_clu_dataset)

#####
#k=4
#####
distance <- dist(new_clu_dataset)
km_4 <- kmeans(new_clu_dataset, 4, nstart = 20)
si_kmeans <- silhouette(km_4$cluster, dist = distance)
summary(si_kmeans)
#silhouette plot
plot(si_kmeans, col=2:5)
tk_4 <- table(km_4$cluster)
#cluster plot
fviz_cluster(km_4, data = new_clu_dataset)

#####
#k=5
#####
distance <- dist(new_clu_dataset)
km_5 <- kmeans(new_clu_dataset, 5, nstart = 20)
si_kmeans <- silhouette(km_5$cluster, dist = distance)
summary(si_kmeans)
#silhouette plot
plot(si_kmeans, col=2:6)
tk_5 <- table(km_5$cluster)
#cluster plot
fviz_cluster(km_5, data = new_clu_dataset)

##Conclusion: 
##From the plots above, it clearly shows that when k = 3 the clustering express best.
##This conclusion can be some different with the PAM50 clinical samples whose k = 4.
##We think this difference might due to the limited samples. Actually we only got 35
##PAM50 genes from 77_cancer_proteomes tables. 

#####add new label to data frame
clu_dataset_patients <- clu_dataset
clu_dataset_patients$Lable4 <- km_4$cluster #k=3 best, but we set k=4 so as to compare with original data
clu_dataset_patients$Lable3 <- km_3$cluster
write.csv(clu_dataset_patients, 'cluster_k_3_4.csv')
######merge by 'Complete TCGA ID'
df3 <- clu_dataset_patients
df4 <- read.csv('clinical_data_breast_cancer.csv')
#rename df4
df4 <- df4[,-1]
df4 <- df4[,-31]
colnames(df4)[1] <- 'Complete TCGA ID'
names(df4)
##merge data
new_clinical <- merge(df4, df3, by = 'Complete TCGA ID')
write.csv(new_clinical, 'new_clinical.csv')
#comparing
#result <- table(new_clinical$Integrated.Clusters..with.PAM50., new_clinical$Lable)
#confusionMatrix(result)

#plot_original <- new_clinical[,c(1, 21, 66)]
#plot()


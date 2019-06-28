library(pathClass)
library(igraph)
setwd("/data/Frank_Huang_lab/net_svm_testing/NBSBM/packedCodes/")
# source('CrossValidation_RFEE.R')
# prostate_cancer_test <- read.table("my_gdsc_dasa_cell_line.txt",check.names=F,header=T)
# labels_test <- as.matrix(read.delim("sensi_cell_resis.txt",header=F))[,2]
# labels_test <- as.factor(labels_test)
# labels_test <- labels_test[-1]

source("nbsbc_16.R")
source("nbsbmnet_revise_LH.R")

prostate_cancer_test <- read.table("./data/my_gdsc_dasa_cell_line.txt",check.names=F,header=T)
labels_test <- as.matrix(read.delim("./data/sensi_cell_resis.txt",header=F))[,2]
labels_test <- labels_test[-1]
labels_test <- as.numeric(labels_test)
labels_test <- as.factor(labels_test)
prostate_cancer_test <- as.matrix(prostate_cancer_test)
network <- as.matrix(read.table("./data/network_ppi_expr_prostate_compare.txt"))
match_index_one <- as.matrix(match(network[,1],rownames(prostate_cancer_test)))
match_index_two <- as.matrix(match(network[,3],rownames(prostate_cancer_test)))
network[,1] <- match_index_one
network[,3] <- match_index_two

fold = 5
repeats = 5
nbsbm.res <- crossval.nbsbm(t(prostate_cancer_test), labels_test, folds=fold, repeats=repeats, parallel=T)
filename <-paste("Nbsbm_whole_cell_line_",fold,"fold",repeats,"repeat_7grid_a1b__1")
plot(nbsbm.res,toFile=T,fname=filename)  


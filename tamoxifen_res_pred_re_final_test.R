setwd("/data/Frank_Huang_lab/net_svm_testing/NBSBM/packedCodes/")

library(pathClass)
library(igraph)
source("nbsbc.R")
source("nbsbmnet_revise_LH.R")

# breast_tamoxifen <- read.delim("breast_tamoxifen_res_halfsample.txt",sep="\t",check.names=F,header=T)
breast_tamoxifen <- read.delim("./data/breast_tamoxifen_res.txt",sep="\t",check.names=F,header=T)
breast_survival <- read.table("./data/breast_survival.txt",sep="\t")
breast_survival <- t(breast_survival)
tlabels <- as.factor(breast_survival)

breast_tamoxifen <- as.matrix(breast_tamoxifen)
network <- as.matrix(read.table("./data/network_ppi_expr_prostate_compare.txt"))
match_index_one <- as.matrix(match(network[,1],rownames(breast_tamoxifen)))
match_index_two <- as.matrix(match(network[,3],rownames(breast_tamoxifen)))
network[,1] <- match_index_one
network[,3] <- match_index_two

fold = 5
repeats = 2
a = 5
b = 9
len = 7

nbsbm.res <- crossval.nbsbm(t(breast_tamoxifen), tlabels, folds=fold, repeats=repeats, parallel=T, a0 = a, b0 = b, len = len)
filename <-paste("testsing_Nbsbm_103",fold,"fold",repeats,"repeat",len,"grid_a",a,"b",b,"0",sep = "_")
plot(nbsbm.res,toFile=T,fname=filename) 


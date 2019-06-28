# source("networkBasedSVM.R")
# source("RRFE.R")
# source('GeneRank.R')
# source('CrossValidation_RFEE.R')
# source("RecursiveFeatureElimination.R")
# source("PathwayMethods.R")
# source("GraphSVM.R")
# source("SVMs.R")
# source("SpanBound.R")
library(pathClass)
library(igraph)

setwd("/data/Frank_Huang_lab/net_svm_testing/NBSBM/packedCodes/")
prostate_res <- read.delim("./data/my_prostate_cell_line_16_median.txt",header=T,check.names=F)
prostate_cell <- c("PC3","DU145","LNCaP","22Rv","WPMY1","VCaP","MDAPCa2b","HPV7","HPV10","RWPE1","RWPE2","NB11","W99","PWR1E","DUCaP","NB26")
class_label <- c(1,1,1,-1,-1,-1,-1,1,1,1,1,1,1,1,-1,1)

protein_inter <- as.matrix(read.table("./data/network_ppi_expr_prostate_compare.txt"))
sub_protein_net <- as.matrix(cbind(protein_inter[,1],protein_inter[,3]))
m <- apply(prostate_res, 2, mean)
s <- apply(prostate_res, 2, sd)
scaledTrain <- (prostate_res - matrix(m, nrow = nrow(prostate_res), ncol = ncol(prostate_res), byrow = T)) /
			matrix(s, nrow = nrow(prostate_res), ncol = ncol(prostate_res), byrow = T)
bio.net <- as.matrix(sub_protein_net)
protein.graph <- graph.data.frame(unique(bio.net))
ad.pro.matrix <- get.adjacency(protein.graph,sparse=F)
get.map <- rownames(prostate_res)
get.map <- unique(get.map)

pro.map <- cbind(get.map,get.map)
colnames(pro.map) <- c("probesetID","graphID")
pro.map <- as.data.frame(pro.map)

source("nbsbc_16.R")
source("nbsbmnet_revise_LH.R")

prostate_cancer_train <- read.table("./data/my_prostate_cell_line_16_median.txt",check.names=F,header=T)
# prostate_cancer_train <- read.table("my_prostate_cell_16.txt",check.names=F,header=T)
labels_train <- c(1,1,1,-1,-1,-1,-1,1,1,1,1,1,1,1,-1,1)
labels_train <- as.factor(labels_train)
network <- as.matrix(read.table("./data/network_ppi_expr_prostate_compare.txt"))
match_index_one <- as.matrix(match(network[,1],rownames(prostate_cancer_train)))
match_index_two <- as.matrix(match(network[,3],rownames(prostate_cancer_train)))
network[,1] <- match_index_one
network[,3] <- match_index_two

nbsbm.res <- crossval.nbsbm(t(prostate_cancer_train), labels_train, folds=3, repeats=5, parallel=T)
filename <-paste("Nbsbm_16_3forld_5repeat_7grid_a1b5_1")
plot(nbsbm.res,toFile=T,fname=filename)




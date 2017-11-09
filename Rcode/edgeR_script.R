#########################edgeR DE##########################################################################
library(limma)
library(edgeR)
library(ggplot2)
library(data.table)

setwd("~/Documents/Young Data/order/")

read.csv("DG.csv", row.names = 1) -> DG
read.csv("SG.csv", row.names = 1) -> SG
read.csv("PC.csv", row.names = 1) -> PC
read.csv("NC.csv", row.names = 1) -> NC

########### DG vs SG ########################################################################################

cbind(DG, SG) -> X

#create edgeR object
group <- rep(c(0, 1), each=3)
Y<- DGEList(counts = X, group = group)

#Filter for lowly expressed genes
keep <- rowSums(cpm(Y) >20) > 3
Y <- Y[keep,]

#give colnames
Y$samples$lib.size <- colSums(Y$counts)

#normalisation
Y <- calcNormFactors(Y, method = "TMM")

#QC
plotMDS(Y)

#design matrix
des <- model.matrix(~1+group)
rownames(des) <- colnames(Y)

#estimate DE
Y <- estimateGLMCommonDisp(Y)
#Y <- estimateGLMTrendedDisp(Y)
#Y <- estimateGLMTagwiseDisp(Y)

#GLM model
fit <- glmFit(Y, des)
lrt <- glmLRT(fit)

#results
setwd("~/Documents/Young Data/DE/")
topTags(lrt, n=Inf) -> DG_SG
write.csv(DG_SG, "DG_SG.txt")

################DG vs PC####################################################################################
cbind(DG, PC) -> X

#create edgeR object
group <- rep(c(0, 1), each=3)
Y<- DGEList(counts = X, group = group)

#Filter for lowly expressed genes
keep <- rowSums(cpm(Y) >20) > 3
Y <- Y[keep,]

#give colnames
Y$samples$lib.size <- colSums(Y$counts)

#normalisation
Y <- calcNormFactors(Y, method = "TMM")

#QC
plotMDS(Y)

#design matrix
des <- model.matrix(~1+group)
rownames(des) <- colnames(Y)

#estimate DE
Y <- estimateGLMCommonDisp(Y)
#Y <- estimateGLMTrendedDisp(Y)
#Y <- estimateGLMTagwiseDisp(Y)

#GLM model
fit <- glmFit(Y, des)
lrt <- glmLRT(fit)

topTags(lrt)
#results
setwd("~/Documents/Young Data/DE/")
topTags(lrt, n=Inf) -> DG_PC
write.csv(DG_PC, "DG_PC.txt")

#############DG vs NC ######################################################################################
cbind(DG, NC) -> X

#create edgeR object
group <- rep(c(0, 1), each=3)
Y<- DGEList(counts = X, group = group)

#Filter for lowly expressed genes
keep <- rowSums(cpm(Y) >20) > 3
Y <- Y[keep,]

#give colnames
Y$samples$lib.size <- colSums(Y$counts)

#normalisation
Y <- calcNormFactors(Y, method = "TMM")

#QC
plotMDS(Y)

#design matrix
des <- model.matrix(~1+group)
rownames(des) <- colnames(Y)

#estimate DE
Y <- estimateGLMCommonDisp(Y)
#Y <- estimateGLMTrendedDisp(Y)
#Y <- estimateGLMTagwiseDisp(Y)

#GLM model
fit <- glmFit(Y, des)
lrt <- glmLRT(fit)

topTags(lrt)
#results
setwd("~/Documents/Young Data/DE/")
topTags(lrt, n=Inf) -> DG_NC
write.csv(DG_NC, "DG_NC.txt")

############# SG vs DG ###################################################################################
cbind(SG, DG) -> X

#create edgeR object
group <- rep(c(0, 1), each=3)
Y<- DGEList(counts = X, group = group)

#Filter for lowly expressed genes
keep <- rowSums(cpm(Y) >20) > 3
Y <- Y[keep,]

#give colnames
Y$samples$lib.size <- colSums(Y$counts)

#normalisation
Y <- calcNormFactors(Y, method = "TMM")

#QC
plotMDS(Y)

#design matrix
des <- model.matrix(~1+group)
rownames(des) <- colnames(Y)

#estimate DE
Y <- estimateGLMCommonDisp(Y)
#Y <- estimateGLMTrendedDisp(Y)
#Y <- estimateGLMTagwiseDisp(Y)

#GLM model
fit <- glmFit(Y, des)
lrt <- glmLRT(fit)

topTags(lrt)
#results
setwd("~/Documents/Young Data/DE/")
topTags(lrt, n=Inf) -> SG_DG
write.csv(SG_DG, "SG_DG.txt")

############ SG vs PC ####################################################################################
cbind(SG, PC) -> X

#create edgeR object
group <- rep(c(0, 1), each=3)
Y<- DGEList(counts = X, group = group)

#Filter for lowly expressed genes
keep <- rowSums(cpm(Y) >20) > 3
Y <- Y[keep,]

#give colnames
Y$samples$lib.size <- colSums(Y$counts)

#normalisation
Y <- calcNormFactors(Y, method = "TMM")

#QC
plotMDS(Y)

#design matrix
des <- model.matrix(~1+group)
rownames(des) <- colnames(Y)

#estimate DE
Y <- estimateGLMCommonDisp(Y)
#Y <- estimateGLMTrendedDisp(Y)
#Y <- estimateGLMTagwiseDisp(Y)

#GLM model
fit <- glmFit(Y, des)
lrt <- glmLRT(fit)

topTags(lrt)
#results
setwd("~/Documents/Young Data/DE/")
topTags(lrt, n=Inf) -> SG_PC
write.csv(SG_PC, "SG_PC.txt")


########### SG vs NC ######################################################################################
cbind(SG, NC) -> X

#create edgeR object
group <- rep(c(0, 1), each=3)
Y<- DGEList(counts = X, group = group)

#Filter for lowly expressed genes
keep <- rowSums(cpm(Y) >20) > 3
Y <- Y[keep,]

#give colnames
Y$samples$lib.size <- colSums(Y$counts)

#normalisation
Y <- calcNormFactors(Y, method = "TMM")

#QC
plotMDS(Y)

#design matrix
des <- model.matrix(~1+group)
rownames(des) <- colnames(Y)

#estimate DE
Y <- estimateGLMCommonDisp(Y)
#Y <- estimateGLMTrendedDisp(Y)
#Y <- estimateGLMTagwiseDisp(Y)

#GLM model
fit <- glmFit(Y, des)
lrt <- glmLRT(fit)

topTags(lrt)
#results
setwd("~/Documents/Young Data/DE/")
topTags(lrt, n=Inf) -> SG_NC
write.csv(SG_NC, "SG_NC.txt")

#############PC  vs  DG ###################################################################################
cbind(PC, DG) -> X

#create edgeR object
group <- rep(c(0, 1), each=3)
Y<- DGEList(counts = X, group = group)

#Filter for lowly expressed genes
keep <- rowSums(cpm(Y) >20) > 3
Y <- Y[keep,]

#give colnames
Y$samples$lib.size <- colSums(Y$counts)

#normalisation
Y <- calcNormFactors(Y, method = "TMM")

#QC
plotMDS(Y)

#design matrix
des <- model.matrix(~1+group)
rownames(des) <- colnames(Y)

#estimate DE
Y <- estimateGLMCommonDisp(Y)
#Y <- estimateGLMTrendedDisp(Y)
#Y <- estimateGLMTagwiseDisp(Y)

#GLM model
fit <- glmFit(Y, des)
lrt <- glmLRT(fit)

topTags(lrt)
#results
setwd("~/Documents/Young Data/DE/")
topTags(lrt, n=Inf) -> PC_DG
write.csv(PC_DG, "PC_DG.txt")

############PC vs SG#####################################################################################
cbind(PC, SG) -> X

#create edgeR object
group <- rep(c(0, 1), each=3)
Y<- DGEList(counts = X, group = group)

#Filter for lowly expressed genes
keep <- rowSums(cpm(Y) >20) > 3
Y <- Y[keep,]

#give colnames
Y$samples$lib.size <- colSums(Y$counts)

#normalisation
Y <- calcNormFactors(Y, method = "TMM")

#QC
plotMDS(Y)

#design matrix
des <- model.matrix(~1+group)
rownames(des) <- colnames(Y)

#estimate DE
Y <- estimateGLMCommonDisp(Y)
#Y <- estimateGLMTrendedDisp(Y)
#Y <- estimateGLMTagwiseDisp(Y)

#GLM model
fit <- glmFit(Y, des)
lrt <- glmLRT(fit)

topTags(lrt)
#results
setwd("~/Documents/Young Data/DE/")
topTags(lrt, n=Inf) -> PC_SG
write.csv(PC_SG, "PC_SG.txt")

###########PC vs NC ######################################################################################
cbind(PC, NC) -> X

#create edgeR object
group <- rep(c(0, 1), each=3)
Y<- DGEList(counts = X, group = group)

#Filter for lowly expressed genes
keep <- rowSums(cpm(Y) >20) > 3
Y <- Y[keep,]

#give colnames
Y$samples$lib.size <- colSums(Y$counts)

#normalisation
Y <- calcNormFactors(Y, method = "TMM")

#QC
plotMDS(Y)

#design matrix
des <- model.matrix(~1+group)
rownames(des) <- colnames(Y)

#estimate DE
Y <- estimateGLMCommonDisp(Y)
#Y <- estimateGLMTrendedDisp(Y)
#Y <- estimateGLMTagwiseDisp(Y)

#GLM model
fit <- glmFit(Y, des)
lrt <- glmLRT(fit)

topTags(lrt)
#results
setwd("~/Documents/Young Data/DE/")
topTags(lrt, n=Inf) -> PC_NC
write.csv(PC_NC, "PC_NC.txt")

############NC vs DG #####################################################################################
cbind(NC, DG) -> X

#create edgeR object
group <- rep(c(0, 1), each=3)
Y<- DGEList(counts = X, group = group)

#Filter for lowly expressed genes
keep <- rowSums(cpm(Y) >20) > 3
Y <- Y[keep,]

#give colnames
Y$samples$lib.size <- colSums(Y$counts)

#normalisation
Y <- calcNormFactors(Y, method = "TMM")

#QC
plotMDS(Y)

#design matrix
des <- model.matrix(~1+group)
rownames(des) <- colnames(Y)

#estimate DE
Y <- estimateGLMCommonDisp(Y)
#Y <- estimateGLMTrendedDisp(Y)
#Y <- estimateGLMTagwiseDisp(Y)

#GLM model
fit <- glmFit(Y, des)
lrt <- glmLRT(fit)

topTags(lrt)
#results
setwd("~/Documents/Young Data/DE/")
topTags(lrt, n=Inf) -> NC_DG
write.csv(NC_DG, "NC_DG.txt")

#########NC  vs  SG #######################################################################################
cbind(NC, SG) -> X

#create edgeR object
group <- rep(c(0, 1), each=3)
Y<- DGEList(counts = X, group = group)

#Filter for lowly expressed genes
keep <- rowSums(cpm(Y) >20) > 3
Y <- Y[keep,]

#give colnames
Y$samples$lib.size <- colSums(Y$counts)

#normalisation
Y <- calcNormFactors(Y, method = "TMM")

#QC
plotMDS(Y)

#design matrix
des <- model.matrix(~1+group)
rownames(des) <- colnames(Y)

#estimate DE
Y <- estimateGLMCommonDisp(Y)
#Y <- estimateGLMTrendedDisp(Y)
#Y <- estimateGLMTagwiseDisp(Y)

#GLM model
fit <- glmFit(Y, des)
lrt <- glmLRT(fit)

topTags(lrt)
#results
setwd("~/Documents/Young Data/DE/")
topTags(lrt, n=Inf) -> NC_SG
write.csv(NC_SG, "NC_SG.txt")
#########NC  vs  SG #######################################################################################
cbind(NC, PC) -> X

#create edgeR object
group <- rep(c(0, 1), each=3)
Y<- DGEList(counts = X, group = group)

#Filter for lowly expressed genes
keep <- rowSums(cpm(Y) >20) > 3
Y <- Y[keep,]

#give colnames
Y$samples$lib.size <- colSums(Y$counts)

#normalisation
Y <- calcNormFactors(Y, method = "TMM")

#QC
plotMDS(Y)

#design matrix
des <- model.matrix(~1+group)
rownames(des) <- colnames(Y)

#estimate DE
Y <- estimateGLMCommonDisp(Y)
#Y <- estimateGLMTrendedDisp(Y)
#Y <- estimateGLMTagwiseDisp(Y)

#GLM model
fit <- glmFit(Y, des)
lrt <- glmLRT(fit)

topTags(lrt)
#results
setwd("~/Documents/Young Data/DE/")
topTags(lrt, n=Inf) -> NC_PC
write.csv(NC_PC, "NC_PC.txt")


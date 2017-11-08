library(limma)
library(matrixStats)
library(edgeR) 
library(splines)

setwd("~/Documents/Jacobs Data/exprs/")
read.csv("genes.csv", row.names = 21) -> mydata
mydata$X <- NULL

rowIQRs(mydata) -> IQ
cbind(mydata, IQ) -> Y

keep <- Y[which(Y$IQ >0.3),]
keep$IQ <- NULL

#plotMDS(mydata)
#plotMDS(keep)

File <- colnames(mydata)
Time <- c(rep(0, 4), c(rep(1, 2)), c(rep(3,3)), c(rep(7, 4)), c(rep(12, 3)), c(rep(18, 3)))
Group <- c(rep(c("Control"),4), rep(c("Treat"),15))
as.array(cbind(File, Group, Time)) -> Targets

spline_time <- ns(Time, df=5)
design <- model.matrix(~Group*spline_time)
rownames(design) <- colnames(mydata)

fit <- lmFit(keep, design)
fit2 <- eBayes(fit)

setwd("~/Documents/Jacobs Data/DE_limma/")
topTable(fit2, coef = 2, n=Inf) -> Z
write.csv(Z, "DE_lim.csv")

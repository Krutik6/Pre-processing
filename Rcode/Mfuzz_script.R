library(Mfuzz)
library(Biobase)

setwd("~/Documents/Denninger Data/normalised data/")
read.csv("datadenn.csv", row.names = 17) -> mdata
mdata$X <- NULL

colnames(mdata) <- c("Control1_0", "Control2_0", "Control3_0", "Arth1_0", "Arth2_0", "Arth3_0",
                     "Arth1_7", "Arth2_7","Arth3_7", "Arth1_14", "Arth2_14", "Arth3_14", 
                     "Dec1_14", "Dec2_14", "Dec3_14")


X <- new("ExpressionSet", exprs=as.matrix(mdata))

X.T <- filter.std(X, min.std = 0)
X.S <- standardise(X)

m1 <- mestimate(X.S)
m1
cl <- mfuzz(X.S, c=16, m=1.158317)
mfuzz.plot(X.S, cl, mfrow = c(4,4))


Ov <- overlap(cl)
PT <- overlap.plot(cl, overlap = Ov, thres = 0.05)

cl_10 <- mfuzz(X.S, c=10, m=1.158317)
mfuzz.plot(X.S, cl_10, mfrow = c(3,4))
O3 <- overlap(cl_10)
overlap.plot(cl_10, overlap = O3, P= PT, thres = 0.05)

cl_25 <- mfuzz(X.S, c=25, m=1.158317)
mfuzz.plot(X.S, cl_25, mfrow = c(5,5))
O4 <- overlap(cl_25)
overlap.plot( cl_25, O4, PT, thres = 0.05)

head(cl_25)
cl_25$centers -> C
cl_25$cluster -> Y
setwd("~/Documents/Denninger Data/mfuzz/")
write.csv(Y, "cluster_data.csv")
cl_25$membership -> Z


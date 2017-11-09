library(lumi)
setwd("~/Documents/Gardiner data/raw data/illumina/")

X <- "GPL6887_MouseWG-6_V2_0_R0_11278593_A.bgx"

Y <-  lumiR(X, convertNuID = FALSE)

library(samr)

setwd("~/Documents/Denninger Data/normalised data/")
read.csv("datadenn.csv", row.names = 17) -> mdata
mdata$X <- NULL


rownames(mdata) -> genes

mdata[,c(4:12)] -> X

y <- paste(rep(1, 9), "Time", rep(c(0, 1, 2), each=3), sep = "")
           
start = 1
for (i in start){ 
  y[i]=paste(y[i], "Start", sep = "")
}
for(i in start+8){
  y[i]=paste(y[i], "End", sep = "")
}

data<- list(x=as.matrix(X),y=y, geneIDs = genes)

####S <- SAM(X, y, resp.type = "One class timecourse", nperms = 200,
         time.summary.type = "slope", genenames = genes, logged2 = FALSE,
         fdr.output = 0.1)
###why SAM not work here????? why samr work.............
S2 <-samr(data = data, resp.type = "One class timecourse", nperms = 100)
delta.table <- samr.compute.delta.table(S2, min.foldchange = 0.1, nvals = 200)
siggenes.table <- samr.compute.siggenes.table(S2, del = 0, data, delta.table = delta.table, all.genes = TRUE)

A <- siggenes.table$genes.up
B <- siggenes.table$genes.lo
c <- rbind(A, B)

lo <- c[as.numeric(c[,8])<1,]


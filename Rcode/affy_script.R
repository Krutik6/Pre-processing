library(affy)

library(mgu74av2.db)


setwd("~/Documents/Jacobs Data/raw/")

read.affybatch("D01.CEL", "D02.CEL", "D03.CEL","D04.CEL","D11.CEL", "D12.CEL", "D31.CEL", "D32.CEL",
               "D33.CEL", "D71.CEL", "D72.CEL", "D73.CEL", "D74.CEL", "D121.CELL", "D122.CEL",
               "D123.CEL", "D181.CEL", "D182.CEL", "D183.CEL",
               description = NULL,
               notes = "", compress = FALSE, rm.mask = F, rm.outliers = F, rm.extra = F,
               verbose = T) -> D0


ND0 <- rma(D0)


setwd("~/Documents/Jacobs Data/exprs/")
write.exprs(ND0, "JDays.csv")

Annot <- data.frame(SYMBOL=sapply(contents(mgu74av2SYMBOL), paste, collapse=","), 
                    DESC=sapply(contents(mgu74av2GENENAME), paste, collapse=","),
                    ENTREZID=sapply(contents(mgu74av2ENTREZID), paste, collapse=","),
                    ENSEMBLID=sapply(contents(mgu74av2ENSEMBL), paste, collapse=","))

Symbol <-  sapply(contents(mgu74av2SYMBOL), paste, collapse= ",")

read.csv("JDays.csv", sep = "\t", row.names = 1) -> X
cbind(X, Symbol) -> Z

head(Symbol, n= 20)

Z[which(!duplicated(Z$Symbol)),] -> Z
Z[!Z$Symbol == "NA",] -> named


write.csv(named, "genes.csv")

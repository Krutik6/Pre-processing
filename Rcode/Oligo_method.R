library(oligo)
library(data.table)
library(mogene10sttranscriptcluster.db)

setwd("~/Documents/Denninger Data/raw data/")

ControlA <- read.celfiles("Control_A.CEL")
ControlB <- read.celfiles("Control_B.CEL")
ControlC <- read.celfiles("Control_C.CEL")
EarlyA <- read.celfiles("Day_0_3_A.CEL")
EarlyB <- read.celfiles("Day_0_3_B.CEL")
EarlyC <- read.celfiles("Day_0_3_C.CEL")
MidA <- read.celfiles("Week_1_2_A.CEL")
MidB <- read.celfiles("Week_1_2_B.CEL")
MidC <- read.celfiles("Week_1_2_C.CEL")
LateA <- read.celfiles("Week_3_4_A.CEL")
LateB <- read.celfiles("Week_3_4_B.CEL")
LateC <- read.celfiles("Week_3_4_C.CEL")
DeclineA <- read.celfiles("Declined_A.CEL")
DeclineB <- read.celfiles("Declined_B.CEL")
DeclineC <- read.celfiles("Declined_C.CEL")

#rma normalisation
ControlA <- rma(ControlA)
ControlB <- rma(ControlB)
ControlC <- rma(ControlC)
EarlyA <- rma(EarlyA)
EarlyB <- rma(EarlyB)
EarlyC <- rma(EarlyC)
MidA <- rma(MidA)
MidB <- rma(MidB)
MidC <- rma(MidC)
LateA <- rma(LateA)
LateB <- rma(LateB)
LateC <- rma(LateC)
DeclineA <- rma(DeclineA)
DeclineB <- rma(DeclineB)
DeclineC <- rma(DeclineC)


#long winded way to get normalised data that is accessable
setwd("~/Documents/Denninger Data/normalised data/")
write.exprs(ControlA, "ControlA.txt")
write.exprs(ControlB, "ControlB.txt")
write.exprs(ControlC, "ControlC.txt")
write.exprs(EarlyA, "EarlyA.txt")
write.exprs(EarlyB, "EarlyB.txt")
write.exprs(EarlyC, "EarlyC.txt")
write.exprs(MidA, "MidA.txt")
write.exprs(MidB, "MidB.txt")
write.exprs(MidC, "MidC.txt")
write.exprs(LateA, "LateA.txt")
write.exprs(LateB, "LateB.txt")
write.exprs(LateC, "LateC.txt")
write.exprs(DeclineA, "DeclineA.txt")
write.exprs(DeclineB, "DeclineB.txt")
write.exprs(DeclineC, "DeclineC.txt")

read.table("ControlA.txt") -> ControlA
read.table("ControlB.txt") -> ControlB
read.table("ControlC.txt") -> ControlC
read.table("EarlyA.txt") -> EarlyA
read.table("EarlyB.txt") -> EarlyB
read.table("EarlyC.txt") -> EarlyC
read.table("MidA.txt") -> MidA
read.table("MidB.txt") -> MidB
read.table("MidC.txt") -> MidC
read.table("LateA.txt") -> LateA
read.table("LateB.txt") -> LateB
read.table("LateC.txt") -> LateC
read.table("DeclineA.txt")-> DeclineA
read.table("DeclineB.txt")-> DeclineB
read.table("DeclineC.txt")-> DeclineC

#get gene entries
Annot <- data.frame(SYMBOL=sapply(contents(mogene10sttranscriptclusterSYMBOL), paste, collapse=","), 
                    DESC=sapply(contents(mogene10sttranscriptclusterGENENAME), paste, collapse=","),
                    ENTREZID=sapply(contents(mogene10sttranscriptclusterENTREZID), paste, collapse=","),
                    ENSEMBLID=sapply(contents(mogene10sttranscriptclusterENSEMBL), paste, collapse=","))

#merge data
cbind(ControlA, ControlB, ControlC, EarlyA, EarlyB, EarlyC
      ,MidA, MidB, MidC, LateA, LateB, LateC, DeclineA, DeclineB
      , DeclineC) -> mydata

cbind(mydata, Annot$SYMBOL) -> anotdata
anotdata[!anotdata$`Annot$SYMBOL` == "NA",] -> named

#remove duplicates
single <- named[!duplicated(named[, 16]),]


write.csv(single, "datadenn.csv")

library(limma)
library(edgeR)
library(data.table)

#Get Time course data
setwd("~/Documents/Denninger Data/normalised data/")
read.csv("datadenn.csv", row.names = 17) -> mydata
mydata$X <- NULL
colnames(mydata) <- c("Con", "Con", "Con", "Early","Early", "Early",
                      "Mid", "Mid", "Mid", "Late", "Late", "Late", "Dec"
                      ,"Dec", "Dec")


Time <- c("Con", "Con", "Con", "Early","Early", "Early",
         "Mid", "Mid", "Mid", "Late", "Late", "Late", "Dec"
         ,"Dec", "Dec")

#MDS
mydata <- as.matrix(mydata)
plotMDS(mydata)


#constuct design matrix
n <- 3
group <-  rep(c("Con", "Early", "Mid", "Late", "Dec"), each= n)
design <- model.matrix(~0+group)
rownames(design) <- Time

#perform DE
fit <- lmFit(mydata, design = design, coef=2)


# make contrasts
Con_early <- makeContrasts("groupCon-groupEarly", levels = design)
Con_Mid <- makeContrasts("groupCon-groupMid", levels = design)
Con_Late <- makeContrasts("groupCon-groupLate", levels = design)
Early_Mid <- makeContrasts("groupEarly-groupMid", levels = design)
Mid_late <- makeContrasts("groupMid-groupLate", levels = design)
Late_Dec <- makeContrasts("groupLate-groupDec", levels = design)
Con_Dec <- makeContrasts("groupCon-groupDec", levels = design)

#apply contrasts onto fit
fitCE <- contrasts.fit(fit, Con_early)
fitCM <- contrasts.fit(fit, Con_Mid)
fitCL <- contrasts.fit(fit, Con_Late)
fitEM <- contrasts.fit(fit, Early_Mid)
fitML <- contrasts.fit(fit, Mid_late)
fitLD <- contrasts.fit(fit, Late_Dec)
fitCD <- contrasts.fit(fit, Con_Dec)

#apply eBayes 
ECE <- eBayes(fitCE, proportion =0.01, stdev.coef.lim =c(0.1, 4), 
              trend = F, robust = F)
ECM <- eBayes(fitCM, proportion =0.01, stdev.coef.lim =c(0.1, 4), 
              trend = F, robust = F)
ECL <- eBayes(fitCL, proportion =0.01, stdev.coef.lim =c(0.1, 4), 
              trend = F, robust = F)
EEM <- eBayes(fitEM, proportion =0.01, stdev.coef.lim =c(0.1, 4), 
              trend = F, robust = F)
EML <- eBayes(fitML, proportion =0.01, stdev.coef.lim =c(0.1, 4), 
        trend = F, robust = F)
ELD <- eBayes(fitLD, proportion =0.01, stdev.coef.lim =c(0.1, 4), 
              trend = F, robust = F)
ECD <- eBayes(fitLD, proportion =0.01, stdev.coef.lim =c(0.1, 4), 
             trend = F, robust = F) 

#format data with topTable
topTable(ECE, coef = 1, adjust.method = "BH", n=Inf) -> CE
topTable(ECM, coef = 1, adjust.method = "BH", n=Inf) -> CM
topTable(ECL, coef = 1, adjust.method = "BH", n=Inf) -> CL
topTable(EEM, coef = 1, adjust.method = "BH", n=Inf) -> EM
topTable(EML, coef = 1, adjust.method = "BH", n=Inf) -> ML
topTable(ELD, coef = 1, adjust.method = "BH", n=Inf) -> LD
topTable(ECD, coef = 1, adjust.method = "BH", n=Inf) -> CD


#write data in new directory
setwd("~/Documents/Denninger Data/differential expression/")
write.csv(CE, "Control_Early.csv")
write.csv(CM, "Control_Mid.csv")
write.csv(CL, "Control_Late.csv")
write.csv(EM, "Early_Middle.csv")
write.csv(ML, "Middle_Late.csv")
write.csv(LD, "Late_Decline.csv")
write.csv(CD, "Control_Decline.csv")



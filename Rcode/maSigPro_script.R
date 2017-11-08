library(maSigPro)

setwd("~/Documents/Denninger Data/normalised data/")
read.csv("datadenn.csv", row.names = 17) -> mdata
mdata$X <- NULL

colnames(mdata) <- c("Control_0D_1", "Control_0D_2", "Control_0D_3", "Arth_0D_1", "Arth_0D_2", "Arth_0D_3",
                     "Arth_7D_1", "Arth_7D_2","Arth_7D_3", "Arth_14D_1", "Arth_14D_2", "Arth_14D_3", 
                     "Dec_14D_1", "Dec_14D_2", "Dec_14D_3")


#create matrix
Time <- c(0,0,0,0,0,0,7,7,7,14,14,14,14,14,14)
replicate <- rep(c(1, 2, 3, 4, 5), each =3)
Control <- c(1,1,1,0,0,0,0,0,0,0,0,0,0,0,0)
Arth <- c(0,0,0,1,1,1,1,1,1,1,1,1,0,0,0)
Dec <- c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1)

Time <- rep(c(0, 7, 14), each =3)
replicate <- rep(c(1:3), each =3)
G <- rep(1, 9)


ss_edes <- cbind(Time, replicate, G)
colnames(mdata) -> X
X[4:12] -> Y
rownames(ss_edes) <- Y


mdata[,4:12] -> x

mas <- maSigPro(x, ss_edes, vars = "each")

mas$sig.genes$independ -> A
View(A$sig.profiles)












ss.GENE <- function(n, r, var11 = 0.01, var12 = 0.02, var13 = 0.02,
                    var14 = 0.02, a1 = 0, a2 = 0, a3 = 0, a4 = 0) {
  tc.dat <- NULL
  for (i in 1:n) {
    gene <- c(rnorm(r, a1, var11), rnorm(r, a1, var12),
              rnorm(r, a3, var13), rnorm(r, a4, var14))
    tc.dat <- rbind(tc.dat, gene)
  }
  tc.dat }



flat <- ss.GENE(n =85, r=3)
induc <- ss.GENE(n = 5, r = 3, a1 = 0, a2 = 0.2, a3 = 0.6, a4 = 1) 
sat <- ss.GENE(n = 5, r = 3, a1 = 0, a2 = 1, a3 = 1.1, a4 = 1.2)
ord <- ss.GENE(n = 5, r = 3, a1 = -0.8, a2 = -1, a3 = -1.3, a4 =-0.9)
ss.DATA <- rbind(flat, induc,sat,ord)






cbind(Time, replicate, G) -> edes
rownames(edes) <- colnames(mdata)

#desing matrix
design <- make.design.matrix(edes, degree=2)
design$groups.vector

# find DE genes
fit <- p.vector(mdata, design, MT.adjust = "BH", Q = 0.05)
fit$i
tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)

#get data
setwd("~/Documents/Denninger Data/maSigPro/")

sigs <-get.siggenes(tstep, rsq = 0.6, vars = "each")
sigs$summary -> frams
write.csv(frams, "list.DE.csv")

sigs$sig.genes$Control -> cont 
View(cont$sig.profiles)

#plots
suma2Venn(sigs$summary[, c(1:3)])

see.genes(sigs$sig.genes$ArthvsControl, show.fit = T, dis = design$dis,
          cluster.method = "hclust", cluster.data = 1, k.mclust = TRUE)

#####try munti time point method
library(splines)

Time <- c("1", "1", "1", "3","3", "3",
          "14", "14", "14", "21", "21", "21", "21"
          ,"21", "21")

Filenames <- c("Con1", "Con2", "Con3", "Early1","Early2", "Early3",
               "Mid1", "Mid2", "Mid3", "Late1", "Late2", "Late3", "Dec1",
               "Dec2", "Dec3")

Group <- c("Control", "Control", "Control", "Disease", "Disease", "Disease",
           "Disease", "Disease", "Disease", "Disease", "Disease"
           , "Disease", "Control", "Control", "Control")

cbind(Filenames, Group, Time) -> X
X <- as.data.table(X)
rownames(X) <- X$Filenames 
X$Filenames <- NULL
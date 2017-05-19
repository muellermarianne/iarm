path <- "H:/R/"
library(foreign)
amts <- read.spss(paste(path,"Datasets/amts/amts.sav",sep=""),to.data.frame=TRUE, use.missings=TRUE)
save(amts,file="amts.Rdata")

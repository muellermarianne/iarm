library(foreign)
desc2 <- read.spss("H:/ERTG/Karlstad/Desc2/DESC2.sav",to.data.frame=TRUE, use.missings=TRUE)
for (i in 5:14) desc2[,i] <- as.numeric(desc2[,i]) -1  # items must be numerical vectors, range 0-4
save(desc2,file="desc2.rdata")

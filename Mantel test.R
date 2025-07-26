library(poppr)
library(adegenet)
setwd("C:/Users/willh/Documents/Thesis Related/DAPC/DAPC Files GU/DAPC")

pop2<- read.fstat("Bobcat_2pop.txt.DAT") #original microsatellite data with all bobcats even without location info
anyNA(pop2@tab)
which(!complete.cases(pop2@tab))

#nothing in pop2 called indNames, won't add names not sure if this is needed
# names <- read.csv("samp_names.csv")
# pop2$indNames <- names



coord <- read.csv("nei_coord2.csv")

pop2$other$xy <- coord  #import coordinates

pop3 <- missingno(pop2, type = "genotype", cutoff = 0)

Dgeo <- dist(pop3$other$xy)
nei1 <- nei.dist(pop3)

#class(nei1)
#class(Dgeo)

ibd <- mantel.randtest(nei1,Dgeo)
ibd

plot(r1 <- mantel.randtest(nei1,Dgeo), main = "Mantel's test")

?mantel.randtest

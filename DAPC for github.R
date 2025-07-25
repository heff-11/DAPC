
#------------------------------------------------------
# Discriminant Analysis of Principal Components (DAPC)
# & Genetic Analyses for Population Structure
#------------------------------------------------------

# Load Packages
install.packages("adegenet")
library(adegenet)

# Set working directory
setwd("C:/Users/willh/Documents/Thesis Related/DAPC/DAPC Files GU/DAPC")


#pop <- read.csv("DAPC Alleles only.csv")

# Load genotype data from .DAT file
pop2<- read.fstat("Bobcat_2pop.txt.DAT")
head(pop)
summary(pop)
#head(pop)
head(pop2)

#data1 <- df2genind(pop, sep= "/,", NA.char= "-9"  ) #turning excel file into "genind" file for adegenet

#is.genind(data1) #confirming the file converted
is.genind(pop2)

#head(data1)

#grp <- find.clusters(data1, max.n.clust = 5)  #This shows PC VAR
#75  #asks how many to retain, then spits out new graph
#2  #this is choosing the number of clusters to go with

#90
#2
#grp


# Perform K-means clustering to find clusters
grp2<-find.clusters(pop2, max.n.clust = 5)
45
3
grp2

#Now run the output from this through DAPC

#dapc1 <- dapc(data1, grp$grp )
#90
#1

#scatter(dapc1)
#compoplot(dapc1, posi='bottomright', lab='', xlab="individuals")

# Run DAPC on clustered data
dapc2<- dapc(pop2, grp2$grp)
45
2

# Plot results
scatter(dapc2)
compoplot(dapc2, posi='bottomright', lab='', xlab="individuals")

# Trying different cluster criteria
grp3<- find.clusters(pop2, criterion = c("diffNgroup"), max.n.clust = 10)
40
4
grp3
dapc3<- dapc(pop2, grp3$grp)
35
4
scatter(dapc3)
compoplot(dapc3, posi='bottomright', lab='', xlab="individuals")
dapc3


#----------------------------------------
#Using optim.a.score to choose the # of PCs
#----------------------------------------

?optim.a.score
optim.a.score(dapc2) # suggests optimal number of PCs=7
dapc2.2<- dapc(pop2, max.n.clust=5)
12
1

scatter(dapc2.2)
comp2.2 <- compoplot(dapc2.2, posi='bottomright', lab='', xlab="individuals")
comp2.2
dapc2.2
write.csv(dapc2.2$posterior, "dapc2.2.csv")  #creating csv w/ probabilistic assignments
dapc2.2$assign

# Repeat for DAPC grp3
optim.a.score(dapc3) #optimal number of PCs=7
dapc3.2<- dapc(pop2, max.n.clust=5)
7
1
scatter(dapc3.2)
compoplot(dapc3.2, posi='bottomright', lab='', xlab="individuals")
write.csv(dapc3.2$posterior, "dapc3.2.csv") #creating csv w/ probabilistic assignments


dapc2.3<-dapc(pop2, max.n.clust=5)
40
1
optim.a.score(dapc2.3) #optimal number of PCs=1
dapc2.4 <- dapc(pop2,max.n.clust=5)
1 # looks like 40 from graph, why is optim.a.score saying 1?
1
scatter(dapc2.4)
compoplot(dapc2.4, posi='bottomright', lab='', xlab="individuals")


citation()


######################################
# Nei's Distance individual
####################################


install.packages('poppr')
library(poppr)

# Reload data for distance analysis
pop2<- read.fstat("Bobcat_2pop.txt.DAT")
pop2

nei.dist(pop2.3, warning = TRUE)


# Tab-delimited genotype data
pop2.3 <- read.delim("Bobcat_2_2pop.txt")

# Importing sample names
names <- read.csv("samp_names.csv")
pop2.3$indNames <- names

# Filtering out individuals with known data issues
pop2.1 <- pop2.3[indNames(pop2) != c("sampleId39", "sampleId38", "sampleId37", "sampleId36", "sampleId35",
               "sampleId34", "sampleId33", "sampleId32", "sampleId31",
               "sampleId30", "sampleId29", "sampleId26")]

pop3 <- read.delim("Bobcat_3pop.DAT") #new file to remove unknown coordinates

coord <- read.csv("nei_coord.csv")

# Compute Nei's genetic distance
nei1 <- nei.dist(pop2, warning = TRUE)

# Assign spatial coordinates
pop3$other$xy <- coord  #import coordinates

# Calculate Euclidean Distance Matrix
Dgeo <- dist(pop3$other$xy) 

#class(nei1)
#class(Dgeo)

# Mantel Test to compare genetic vs. geographic distance (testing isolation by distance)
ibd <- mantel.randtest(nei1,Dgeo) 
ibd


pop2.1 <- as.matrix(pop2)

?mantel.randtest


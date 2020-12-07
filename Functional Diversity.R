install.packages("FD")
library(FD)

setwd("~/Dropbox/Paper II - Seagrass diversity fisheries/Data/FD")

seagrass.trait<- read.csv("species_trait_matrix.csv", row.names=1, sep=";")
seagrass.abun<- as.matrix(read.csv("species_abundance_matrix.csv", row.names=1, sep=";"))
seagrass.trait
seagrass.abun

dummy$abun

# examples

ex1 <- gowdis(seagrass.trait)
ex1

ex2 <- functcomp(seagrass.trait, seagrass.abun)
ex2

ex3 <- dbFD(seagrass.trait, seagrass.abun)
ex3
data.frame(ex3)
write.csv(ex3, "FD.csv")

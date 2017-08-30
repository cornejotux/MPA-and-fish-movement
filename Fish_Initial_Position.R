### Author: Jorge Cornejo
### Date: March 13, 2013
### Goal: Define position of the fish for the initial
#     time step of the model.

setwd("~/Model/")
data <- read.table(file="2013initialposwithangles.dat", sep='', dec='.')

## How many Fish??
HowMany <- 60000

data$V1 <- as.character(data$V1)
data$V1[1] <- paste(" ", data$V1[1])
data$V1[2:length(data$V1)] <- paste(" ", data$V1[2:length(data$V1)])

data$V2 <- as.character(data$V2 - 15)
data$V2 <- paste("  ", data$V2)

data$V3 <- as.character(paste("   ", data$V3))



write.table(data, file="FishInitial.dat", sep='', dec='.', row.names=F, col.names=F, quote=F)


x <- runif(n=HowMany, min=0.00001, max=50)
x <- ifelse(x < 0, x*-1, x)
y <- runif(n=HowMany, min=0.00001, max=100)
y <- ifelse(y < 0, y*-1, y)
a <- runif(n=HowMany, min=0.00001, max=359.99999)
data <- data.frame(V1 =x, V2 = y, V3 = a)

data$V1 <- as.character(data$V1)
data$V1[1] <- paste(" ", data$V1[1])
data$V1[2:length(data$V1)] <- paste(" ", data$V1[2:length(data$V1)])

data$V2 <- as.character(data$V2)
data$V2 <- paste("  ", data$V2)

data$V3 <- as.character(paste("   ", data$V3))

write.table(data, file="FishInitialMedium.dat", sep='', dec='.', row.names=F, col.names=F, quote=F)


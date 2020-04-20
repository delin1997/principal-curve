library(R.matlab)
library(rgl)

source("pgh.R")
source("principal_curve_KDE.R")

volume <- readMat('sample.mat')$vol
label  <- readMat('true_label.mat')$md
label[label!=0] = 1

data <- which(label!=0, arr.ind = TRUE)
plot3d(data,col=rgb(0,0,1,0.2))

kernel_sigma <- 1.5
pc_projection <- principal_curve_KDE(data, data, kernel_sigma,1)
plot3d(pc_projection,col=rgb(0,0,1,0.2))
source("pgh.R")
source("principal_curve_KDE.R")
library(rgl)
library(tcltk)

t = seq(0,10*pi,pi/50)
x = 5*sin(t)
y = 5*cos(t)
z = t

x_noise <- x+rnorm(length(t),0,0.5)
y_noise <- y+rnorm(length(t),0,0.5)
z_noise <- z+rnorm(length(t),0,0.5)
data <- cbind(x_noise,y_noise,z_noise)

kernel_sigma <- 1.5

pc_projection <- principal_curve_KDE(data, data, kernel_sigma,1)

plot3d(data,xlab="x", ylab = "y", zlab = "z")
lines3d(pc_projection,col=rgb(0,0,1,0.5))

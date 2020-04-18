principal_curve_KDE <- function(X, Xinit, kernel_sigma, targetdim){
  point <- Xinit
  N <- nrow(X)
  dim <- ncol(X)
  threshold <- 1e-3
  pc <- matrix(0, nrow = N, ncol = dim) 
  for (ind in 1:nrow(point)) {
    flag <- 0
    pghlist <- pgh(point[ind, ], X, N, dim, kernel_sigma)
    point_pc1 <- point[ind, ]
    eigenlist <- eigen(pghlist$SI)
    ConstrainedSpace <- eigenlist$vectors[, 1:(dim-targetdim)]
    gra <- pghlist$g
    H <- pghlist$H
    if(abs(t(gra)%*%H%*%gra/(norm(t(gra)%*%H,"2")*norm(t(gra)%*%H,"2")))<0.01){
      flag <- 1
      pc[ind, ]= point_pc1
    }
    if(!flag){
      for (a in 1:20) {
        G <- kernel_matrix(point_pc1, X, N, dim, kernel_sigma)
        num1 <- rowSums(matrix(rep(G, dim), nrow = dim, byrow = T)*t(X))
        den1 <- sum(G)
        pghlist <- pgh(point_pc1, X, N, dim, kernel_sigma)
        point_pc1_old <- point_pc1
        for (c in 1:ncol(ConstrainedSpace)) {
          direction <- ConstrainedSpace[, c]
          point_pc1 <- point_pc1 + direction%*%(t(direction)%*%(num1/den1-point_pc1))
        }
        if(sum(abs(point_pc1_old-point_pc1)<threshold)==dim){
          break
        }
      }
      pc[ind, ]= point_pc1
    }
  }
  return(pc)
}
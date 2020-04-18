pgh <- function(x, data, N, dim, kernel_sigma){
  Kvector <- kernel_matrix(x, data, N, dim, kernel_sigma)
  Xvector <- matrix(rep(x,N),nrow = dim)-t(data)
  GX <- Xvector*matrix(rep(Kvector,dim), nrow = dim, byrow = T)
  p <- mean(Kvector)
  g <- -rowMeans(GX)/kernel_sigma^2
  H <- 1/(N*(kernel_sigma^4)) * (GX%*%t(Xvector))-diag(dim)*p/(kernel_sigma^2)
  SI <- -(1/p)*H+(1/p)^2*g%*%t(g)
  return(list(p=p, g=g, H=H, SI=SI))
}

kernel_matrix <- function(x, data, N, dim, kernel_sigma){
  C <- 1/((2*pi)^(dim/2)*kernel_sigma^dim)
  dx <- (matrix(rep(x,N),nrow = dim)-t(data))/matrix(rep(kernel_sigma,dim*N),nrow = dim, ncol = N)
  K <- C*exp(-0.5*colSums(dx^2))
  return(K)
}

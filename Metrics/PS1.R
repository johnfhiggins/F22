N = 500

X_sample <- rnorm(N)
eps <- rnorm(N)
theta <- c(0, 1,1)
Y_sample <- as.numeric(X_sample + eps > 0)

max_like(X, theta, N, S)
  Y_sample = 

    
SLL = function(X, Y, theta, z_seq){
  val <- 0.0
  for (i in 1:length(X)){
    p <- p_S(X[i], theta, z_seq)
    val <- val + log(p^(Y[i])(1-p)^(1-Y[i]))/length(X)
  }
  return(val)
}


p_S = function(xi, theta, z_seq){
  val <- 0.0
  for (s in 1:length(z_seq)){
    val <- val + pnorm(theta[1] + (sqrt(theta[3])*z_seq[s] + theta[2])*xi)/length(z_seq)
  }
  return(val)
}
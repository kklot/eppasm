# generate relative partner change rate from uniform
sexual_rate <- function(x, fp){
  # x[c(2, 4, 6, 8)] <- exp(x[c(2, 4, 6, 8)]) # sd and kurtosis are postive
  # cbind(gamlss.dist::dSHASHo(1:66/10, x[1], x[2], x[3], x[4]), 
  #       gamlss.dist::dSHASHo(1:66/10, x[5], x[6], x[7], x[8]))+1e-6
	matrix(x, nrow = 66, ncol = 2)
}

sexual_rate_prior <- function(x, fp)
  sum(dunif(x, 0, 1, log=TRUE))

sexual_rate_sample <- function(n, fp){
  mat <- matrix(NA, n, 66 * 2)
  mat[,] <- runif(66*2*n, 0, 1)
  mat
}

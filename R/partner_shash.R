# generate partner change rate  66x2 from SHASHo
sexual_rate <- function(x, fp){
  #   1:4 male mean sd skew kurtosis
  #   2:8 female ...
  x[c(2, 4, 6, 8)] <- exp(x[c(2, 4, 6, 8)]) # sd and kurtosis are postive
  cbind(gamlss.dist::dSHASHo(1:66/10, x[1], x[2], x[3], x[4]), 
        gamlss.dist::dSHASHo(1:66/10, x[5], x[6], x[7], x[8]))+1e-6
}

sexual_rate_prior <- function(x, fp)
  sum(dnorm(x[c(1,5)], 1, .3, log=TRUE)) +
    sum(dnorm(x[c(2,6)], -1, .5, log=TRUE)) +
      sum(dnorm(x[c(3,7)], 0, .5, log=TRUE)) +
        sum(dnorm(x[c(4,8)], 0, .5, log=TRUE))

sexual_rate_sample <- function(n, fp){
  mat <- matrix(NA, n, 8)
  mat[, c(1,5)] <- rnorm(n*2, 1, .3)
  mat[, c(2,6)] <- rnorm(n*2, -1, .5)
  mat[, c(3,7)] <- rnorm(n*2, 0, .5)
  mat[, c(4,8)] <- rnorm(n*2, 0, .5)
  mat
}

# interpolator
idw <- function(x, D, y){
  d <- (D - x)^2
  if (min(d) == 0) return (y[which.min(d)])
  fhat <- sum(y/d) / sum(1/d)
  return (fhat)
}

# compute c.d.f. for generalized chisq
# by importance sampling
pgchisq <- function(q, lambda, delta, D=NULL){
  p <- length(lambda)
  if (is.null(D)) D <- sp(1000,p,dist.str=rep("normal",p))$sp
  # compute importance weight
  lambda <- p / q * lambda # p to q ratio adjusted lambda
  D.wts <- apply(D, 1, function(x) -sum(dnorm(x,log=T)))
  for (j in 1:p){
    sd <- sqrt(lambda[j])
    mean <- sqrt(delta[j]) * sd
    D.wts <- D.wts + dnorm(D[,j], mean=mean, sd=sd, log=T)
  }
  D.wts <- exp(D.wts)
  D.sum <- apply(D, 1, function(x) sum(x^2))
  prob <- mean(D.wts * (D.sum < p))
  # if (prob < 1e-9) prob <- 0 # truncate at 1e-9
  return (prob)
}


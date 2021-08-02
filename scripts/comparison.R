#########################################################
# scalar output 
rm(list=ls())
set.seed(20210416)
# one dimensional case
f <- function(x){
  v <- sin(30*(x-.9)^4)*cos(2*(x-.9))+(x-.9)/2
  return (v)
}

# goal to find x such that f(x) = 0
target <- 0.15
# initial 5 point designs
n <- 5
D <- 0.5 + 0.5 * cos((2*(n:1)-1)/(2*n)*pi)
y <- f(D)
e <- (y - target)^2
# discretize point for visualization
xs <- seq(0, 1, length.out = 101)
xs <- xs[-51] # remove 0.5 in D
ys <- f(xs)
es <- (ys-target)^2

# visualization
plot(xs, ys, "l", xlab = expression(theta), ylab = expression(H(theta)), lwd = 2)
abline(h = target, lty = 2, col = "blue", lwd = 2)
points(D, y, col = "green", pch = 16, cex = 2)

plot(xs, es, "l", xlab = expression(theta), ylab = expression(f(theta)), lwd = 2)
abline(h = 0, lty = 2, col = "blue", lwd = 2)
points(D, e, col = "green", pch = 16, cex = 2)

# surrogate model by gp
library(mlegp)

#########################################################
# direct ei
# normal distribution approach
# surrogate model
set.seed(20210416)
gp.model <- mlegp(matrix(D,ncol=1),e)
es.pred <- predict.gp(gp.model,matrix(xs,ncol=1),se.fit=TRUE)
es.hat <- es.pred$fit
es.se <- es.pred$se.fit

# compute ei
emin <- min(e)
ei <- rep(0, length(xs))
for (i in 1:length(xs)){
  u <- (emin - es.hat[i]) / es.se[i]
  ei[i] <- es.se[i] * (u * pnorm(u) + dnorm(u))
}
# find the best ei points
D.new <- xs[which.max(ei)]

# ei visualization
plot(xs, ei, "l", xlab = expression(theta), ylab = "EI", lwd = 2)
points(D, rep(0,length(D)), col = "green", pch = 16, cex = 2)
points(D.new, 0, col = "red", pch = 15, cex = 2)

# visualize in e
es.lb <- es.hat + qnorm(0.1) * es.se
es.ub <- es.hat + qnorm(0.9) * es.se
plot(xs, es, "l", xlab = expression(theta), ylab = expression(f(theta)), lwd = 2)
abline(h=0, lty = 2, col = "blue", lwd = 2)
points(D, e, col = "green", pch = 16, cex = 2)
lines(xs, es.hat, col = "red", lty = 2, lwd = 2)
lines(xs, es.lb, col = "red", lty = 3, lwd = 2)
lines(xs, es.ub, col = "red", lty = 3, lwd = 2)
points(D.new, (f(D.new)-target)^2, col = "red", pch = 15, cex = 2)

#########################################################
# log approach
# surrogate model
loge <- log(e)

set.seed(20210416)
gp.model <- mlegp(matrix(D,ncol=1),loge)
loges.pred <- predict.gp(gp.model,matrix(xs,ncol=1),se.fit=TRUE)
loges.hat <- loges.pred$fit
loges.se <- loges.pred$se.fit

# compute ei
emin <- min(e)
logemin <- min(loge)
ei <- rep(0, length(xs))
for (i in 1:length(xs)){
  mu <- loges.hat[i]
  sigma <- loges.se[i]
  ei[i] <- pnorm((logemin-mu)/sigma)*emin - 
    exp(mu+sigma^2/2)*pnorm((logemin-mu-sigma^2)/sigma)
}
# find the best ei points
D.new <- xs[which.max(ei)]

# ei visualization
plot(xs, ei, "l", xlab = expression(theta), ylab = "EI", lwd = 2)
points(D, rep(0,length(D)), col = "green", pch = 16, cex = 2)
points(D.new, 0, col = "red", pch = 15, cex = 2)

# visualize in loge
loges.lb <- loges.hat + qnorm(0.1) * loges.se
loges.ub <- loges.hat + qnorm(0.9) * loges.se
plot(xs, log(es), "l", xlab = expression(theta), ylab = expression(f(theta)), lwd = 2, ylim = c(-8,0))
# abline(h=0, lty = 2, col = "blue", lwd = 2)
points(D, log(e), col = "green", pch = 16, cex = 2)
lines(xs, loges.hat, col = "red", lty = 2, lwd = 2)
lines(xs, loges.lb, col = "red", lty = 3, lwd = 2)
lines(xs, loges.ub, col = "red", lty = 3, lwd = 2)

# visualize in e
es.mean <- exp(loges.hat + loges.se^2/2)
es.lb <- exp(loges.hat + qnorm(0.1) * loges.se)
es.ub <- exp(loges.hat + qnorm(0.9) * loges.se)
plot(xs, es, "l", xlab = expression(theta), ylab = expression(f(theta)), lwd = 2)
abline(h=0, lty = 2, col = "blue", lwd = 2)
points(D, e, col = "green", pch = 16, cex = 2)
lines(xs, es.mean, col = "red", lty = 2, lwd = 2)
lines(xs, es.lb, col = "red", lty = 3, lwd = 2)
lines(xs, es.ub, col = "red", lty = 3, lwd = 2)
points(D.new, (f(D.new)-target)^2, col = "red", pch = 15, cex = 2)

#########################################################
# chi-square approach
# surrogate model
set.seed(20210416)
gp.model <- mlegp(matrix(D,ncol=1),y)
ys.pred <- predict.gp(gp.model,matrix(xs,ncol=1),se.fit=TRUE)
ys.hat <- ys.pred$fit
ys.se <- ys.pred$se.fit

# compute ei
emin <- min(e)
ei <- rep(0, length(xs))
es.hat <- rep(0,length(xs))
es.lb <- rep(0,length(xs))
es.ub <- rep(0,length(xs))
for (i in 1:length(xs)){
  gamma <- ys.se[i]
  delta <- (ys.hat[i]-target)^2/gamma^2
  t <- emin/gamma^2
  K <- 1
  FK <- pchisq(t,df=K,ncp=delta)
  FK2 <- pchisq(t,df=K+2,ncp=delta)
  FK4 <- pchisq(t,df=K+4,ncp=delta)
  ei[i] <- emin*FK - gamma^2*(K*FK2+delta*FK4)
  es.hat[i] <- gamma^2*delta
  es.lb[i] <- gamma^2*qchisq(0.1,df=K,ncp=delta)
  es.ub[i] <- gamma^2*qchisq(0.9,df=K,ncp=delta)
}
# find the best ei points
D.new <- xs[which.max(ei)]

# ei visualization
plot(xs, ei, "l", xlab = expression(theta), ylab = "EI", lwd = 2)
points(D, rep(0,length(D)), col = "green", pch = 16, cex = 2)
points(D.new, 0, col = "red", pch = 15, cex = 2)

# visualize in e
plot(xs, es, "l", xlab = expression(theta), ylab = expression(f(theta)), lwd = 2)
abline(h=0, lty = 2, col = "blue", lwd = 2)
points(D, e, col = "green", pch = 16, cex = 2)
lines(xs, es.hat, col = "red", lty = 2, lwd = 2)
lines(xs, es.lb, col = "red", lty = 3, lwd = 2)
lines(xs, es.ub, col = "red", lty = 3, lwd = 2)
points(D.new, (f(D.new)-target)^2, col = "red", pch = 15, cex = 2)

# visualize in y
plot(xs, ys, "l", lwd = 2,
     xlab = expression(x), ylab = expression(h(x)))
points(D, y, col = "green", pch = 16, cex = 2)
lines(xs, ys.hat, col = "red", lty = 2, lwd = 2)
lines(xs, ys.hat - 2 * ys.se, col = "red", lty = 3, lwd = 2)
lines(xs, ys.hat + 2 * ys.se, col = "red", lty = 3, lwd = 2)
abline(h=target, lty = 2)
points(D.new, f(D.new), col = "red", pch = 16, cex = 1)

#########################################################
#########################################################
#########################################################
# functional output 
rm(list=ls())
set.seed(20210416)
# functional output example
z <- seq(0, 1, length.out = 1001)
f <- function(x){
  zs <- 5 * z
  x1 <- 1 + 2 * x
  x2 <- (1 - x/2)
  ys <- exp(-x2^2/3 * zs / 2) * cos(sqrt(1 + x1^2 - x2^3/4) * zs)
  return (ys)
}

# target data
target.x <- 0.7
target.y <- f(target.x)

# initial 5 point designs
n <- 5
D <- 0.5 + 0.5 * cos((2*(n:1)-1)/(2*n)*pi)
Y <- matrix(NA,nrow=length(z),ncol=length(D))
for (i in 1:length(D)) Y[,i] <- f(D[i])
e <- rep(0,length(D))
for (i in 1:length(D)) e[i] <- mean((Y[,i]-target.y)^2)

# discretize point for visualization
xs <- seq(0, 1, length.out = 101)
xs <- xs[-51] # remove 0.5 in D
Ys <- matrix(NA,nrow=length(z),ncol=length(xs))
for (i in 1:length(xs)) Ys[,i] <- f(xs[i])
matplot(z, Ys, type = "l", xlab = "t", ylab = expression(paste("H(t;",theta,")",sep="")))
lines(z, target.y, lty = 1, lwd = 5)
# compute the mse
es <- rep(0, length(xs))
for (i in 1:length(xs)) es[i] <- mean((Ys[,i]-target.y)^2)

plot(xs, es*length(z), "l", xlab = expression(theta), ylab = expression(f(theta)), lwd = 2)
abline(h=0, lty = 2, col = "blue", lwd = 2)
points(D, e*length(z), pch = 16, cex = 2, col = "green")

#########################################################
# direct ei
# normal distribution approach
# surrogate model
library(mlegp)
set.seed(20210416)
gp.model <- mlegp(matrix(D,ncol=1),e)
es.pred <- predict.gp(gp.model,matrix(xs,ncol=1),se.fit=TRUE)
es.hat <- es.pred$fit
es.se <- es.pred$se.fit

# compute ei
emin <- min(e)
ei <- rep(0, length(xs))
for (i in 1:length(xs)){
  if (es.se[i] > 1e-6){
    u <- (emin - es.hat[i]) / es.se[i]
    ei[i] <- es.se[i] * (u * pnorm(u) + dnorm(u))
  }
}
# find the best ei points
D.new <- xs[which.max(ei)]

# ei visualization
plot(xs, ei, "l", xlab = expression(theta), ylab = "EI", lwd=2)
points(D, rep(0,length(D)), col = "green", pch = 16, cex=2)
points(D.new, 0, col = "red", pch = 15, cex=2)

# visualize in e
# multiply by length(z) gets back to SSE
es.lb <- es.hat + qnorm(0.1) * es.se
es.ub <- es.hat + qnorm(0.9) * es.se
plot(xs, es, "l", xlab = expression(theta), ylab = expression(f(theta)), lwd=2)
points(D, e, col = "green", pch = 16, cex = 2)
lines(xs, es.hat, col = "red", lty = 2, lwd = 2)
lines(xs, es.lb, col = "red", lty = 3, lwd = 2)
lines(xs, es.ub, col = "red", lty = 3, lwd = 2)
points(D.new, mean((f(D.new)-target.y)^2), col = "red", pch = 15, cex = 2)

#########################################################
# log normal distribution approach
# surrogate model
loge <- log(e)
set.seed(20210416)
gp.model <- mlegp(matrix(D,ncol=1),loge)
loges.pred <- predict.gp(gp.model,matrix(xs,ncol=1),se.fit=TRUE)
loges.hat <- loges.pred$fit
loges.se <- loges.pred$se.fit

# compute ei
emin <- min(e)
logemin <- min(loge)
ei <- rep(0, length(xs))
for (i in 1:length(xs)){
  mu <- loges.hat[i]
  sigma <- loges.se[i]
  ei[i] <- pnorm((logemin-mu)/sigma)*emin - 
    exp(mu+sigma^2/2)*pnorm((logemin-mu-sigma^2)/sigma)
}
# find the best ei points
D.new <- xs[which.max(ei)]

# ei visualization
plot(xs, ei, "l", xlab = expression(theta), ylab = "EI", lwd = 2)
points(D, rep(0,length(D)), col = "green", pch = 16, cex = 2)
points(D.new, 0, col = "red", pch = 15, cex = 2)

# visualize in loge
loges.lb <- loges.hat + qnorm(0.1) * loges.se
loges.ub <- loges.hat + qnorm(0.9) * loges.se
plot(xs, log(es), "l", xlab = expression(theta), ylab = expression(f(theta)), lwd = 2, ylim = c(-8,0))
# abline(h=0, lty = 2, col = "blue", lwd = 2)
points(D, log(e), col = "green", pch = 16, cex = 2)
lines(xs, loges.hat, col = "red", lty = 2, lwd = 2)
lines(xs, loges.lb, col = "red", lty = 3, lwd = 2)
lines(xs, loges.ub, col = "red", lty = 3, lwd = 2)

# visualize in e
es.mean <- exp(loges.hat + loges.se^2/2)
es.lb <- exp(loges.hat + qnorm(0.1) * loges.se)
es.ub <- exp(loges.hat + qnorm(0.9) * loges.se)
plot(xs, es, "l", xlab = expression(theta), ylab = expression(f(theta)), lwd = 2)
abline(h=0, lty = 2, col = "blue", lwd = 2)
points(D, e, col = "green", pch = 16, cex = 2)
lines(xs, es.mean, col = "red", lty = 2, lwd = 2)
lines(xs, es.lb, col = "red", lty = 3, lwd = 2)
lines(xs, es.ub, col = "red", lty = 3, lwd = 2)
points(D.new, mean((f(D.new)-target.y)^2), col = "red", pch = 15, cex = 2)

#########################################################
# chi-square approach 
# load library
#setwd("~/gatech/research/project/material/calibration/paper2/scripts/")
source("lib.R") # generalized chi-square
library(fda) # fpca
library(statmod) # quadrature
# library(CompQuadForm)
library(support) # chi-sq ub and lb estimation

# construct fpca by bspline
set.seed(20210416)
nBasis <- 10
basisobj <- create.bspline.basis(range(z),nBasis)
basis <- smooth.basis(z, Y, basisobj)
fpca <- pca.fd(basis$fd, nBasis)
var.pca <- cumsum(fpca$varprop)
nharm <- sum(var.pca < 0.99) + 1
fpca <- pca.fd(basis$fd, nharm)
# obtain eigenvalue
fpca.eig.value <- fpca$values[1:nharm]
# obtain mean and eigenvector
fpca.mean <- eval.fd(z, fpca$meanfd)
fpca.eig.vector <- eval.fd(z, fpca$harmonics)
# obtain standardized pc score (mean 0, var 1)
fpca.std.score <- t(t(fpca$scores)/sqrt(fpca.eig.value))

# compute target standardized pc score
target.basis <- smooth.basis(z,matrix(target.y-fpca.mean,ncol=1),basisobj)
target.score <- inprod(target.basis$fd, fpca$harmonics)
target.std.score <- target.score / sqrt(fpca.eig.value)

# target reconstruction
target.recon <- fpca.mean
for (i in 1:nharm){
  part.recon <- target.std.score[i]*sqrt(fpca.eig.value[i])*fpca.eig.vector[,i]
  target.recon <- target.recon + part.recon
}
mean((target.recon - target.y)^2)
plot(z, target.y, "l")
lines(z, target.recon, col = "red", lty = 2)

# surrogate model
models <- list()
for (i in 1:nharm) models[[i]] <- mlegp(matrix(D,ncol=1), fpca.std.score[,i])

# quadrature points
quad.node <- gauss.quad(10, kind = "legendre")
quad.node$nodes <- (quad.node$nodes+1)/2 # rescale to [0,1]
quad.node$weights <- quad.node$weights/2 # adjust the weight
# support points - chisq approximation
sp <- sp(1000,nharm,dist.str=rep("normal",nharm))$sp

# compute ei
emin <- min(e)
ei <- rep(0, length(xs))
es.hat <- rep(0,length(xs))
es.lb <- rep(0,length(xs))
es.ub <- rep(0,length(xs))
for (j in 1:length(xs)){
  score.fit <- rep(0, nharm)
  score.se <- rep(0, nharm)
  for (i in 1:nharm){
    score <- predict(models[[i]], xs[j], se.fit=TRUE)
    score.fit[i] <- score$fit
    score.se[i] <- score$se.fit
  }
  lambda <- fpca.eig.value * score.se^2 / emin
  delta <- (score.fit - target.std.score)^2 / (score.se)^2
  quad.node.y <- rep(NA, length(quad.node$nodes))
  for (i in 1:length(quad.node$nodes)){
    quad.node.y[i] <- 1 - liu(q=quad.node$nodes[i], lambda=lambda, delta=delta)
  }
  ei[j] <- sum(quad.node.y * quad.node$weights) * emin
  es.hat[j] <- sum(fpca.eig.value * (score.fit - target.std.score)^2)
  es.hat.edf <- rep(0, nrow(sp))
  for (i in 1:nrow(sp)){
    score <- score.fit + score.se * sp[i,]
    es.hat.edf[i] <- sum(fpca.eig.value * (score - target.std.score)^2)
  }
  es.lb[j] <- quantile(es.hat.edf, 0.1)
  es.ub[j] <- quantile(es.hat.edf, 0.9)
}

# find the best ei points
D.new <- xs[which.max(ei)]

# ei visualization
plot(xs, ei, "l", xlab=expression(theta), ylab="EI", lwd = 2)
points(D, rep(0,length(D)), col = "green", pch = 16, cex = 2)
points(D.new, 0, col = "red", pch = 15, cex = 2)

# visualize in e
plot(xs, es, "l", xlab=expression(theta), ylab=expression(f(theta)), lwd = 2)
points(D, e, col = "green", pch = 16, cex = 2)
lines(xs, es.hat, col = "red", lty = 2, lwd = 2)
lines(xs, es.lb, col = "red", lty = 3, lwd = 2)
lines(xs, es.ub, col = "red", lty = 3, lwd = 2)
points(D.new, mean((f(D.new)-target.y)^2), col = "red", pch = 15, cex = 2)

# visualize in y
D.new.score.fit <- rep(0, nharm)
D.new.score.se <- rep(0, nharm)
for (i in 1:nharm){
  score <- gp.predict(D.new, models[[i]])
  D.new.score.fit[i] <- score$fit
  D.new.score.se[i] <- score$se.fit
}
D.new.yhat <- fpca.mean
for (i in 1:nharm){
  part.recon <- D.new.score.fit[i]*sqrt(fpca.eig.value[i])*fpca.eig.vector[,i]
  D.new.yhat <- D.new.yhat + part.recon
}
D.new.yhat.edf <- matrix(NA, nrow = length(z), ncol = nrow(sp))
for (j in 1:nrow(sp)){
  score <- D.new.score.fit + D.new.score.se * sp[j,]
  yhat <- fpca.mean
  for (i in 1:nharm){
    part.recon <- score[i]*sqrt(fpca.eig.value[i])*fpca.eig.vector[,i]
    yhat <- yhat + part.recon
  }
  D.new.yhat.edf[,j] <- yhat 
}
D.new.yhat.lb <- apply(D.new.yhat.edf, 1, function(x) quantile(x,0.1))
D.new.yhat.ub <- apply(D.new.yhat.edf, 1, function(x) quantile(x,0.9))

plot(z, target.y, "l")
lines(z, D.new.yhat, col = "red", lty = 2)
lines(z, D.new.yhat.lb, col = "red", lty = 3)
lines(z, D.new.yhat.ub, col = "red", lty = 3)
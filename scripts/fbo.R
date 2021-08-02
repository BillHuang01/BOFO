rm(list=ls())

# load library
source("lib.R") # generalized chi-square
library(fda) # fpca
library(mlegp) # gaussian process
library(randtoolbox) # sobol point
library(statmod) # quadrature
library(support) # support point
# library(CompQuadForm) # generalized chi-square

# load design
ini <- read.csv("../results/maxpro/maxpro_design.csv")
D <- as.matrix(ini[1:nrow(ini),c(2:6)])
colnames(D) <- NULL
# transform it back to 0 and 1 scale
# Diffusivity, (1e-12,1e-9), log
D[,1] <- (log(D[,1]) - log(1e-12)) / (log(1e-9) - log(1e-12))
# Surface Concentration, (4e-3,5e-3), original
D[,2] <- (D[,2] - 4e-3) / (5e-3 - 4e-3)
# Accessible Polymer Concentration: (5e-3,6e-3), original
D[,3] <- (D[,3] - 5e-3) / (6e-3 - 5e-3)
# Hindering, (500,2500), original
D[,4] <- (D[,4] - 500) / (2500 - 500)
# Reaction rate, (1e-3,1e1), log
D[,5] <- (log(D[,5]) - log(1e-3)) / (log(1e1) - log(1e-3))
summary(D)
pairs(D)

# load experimental data
temp <- 130
torr <- 8.7
length <- 483
setting <- sprintf("%sc_%storr_%snm", temp, torr, length)
sorption <- read.csv(paste("../data/exp_sorption_",setting,".csv",sep=""))
sorption <- sorption[1:125001,] # only consider up to time 62500
desorption <- read.csv(paste("../data/exp_desorption_",setting,".csv",sep=""))
desorption <- desorption[2:115001,] # only consider from 62500 to 120000
exp.t <- seq(0,120000,by=0.5)
exp.y <- c(sorption$mass_uptake,desorption$mass_uptake)
plot(exp.t, exp.y, "l")

# load pde solver data
pde.t <- read.table("../results/timeIndex.txt")[,1]
pde.y <- matrix(NA, nrow = length(pde.t), ncol = nrow(D))
for (i in 1:nrow(D)){
  pde.y[,i] <- read.table(sprintf("../results/maxpro/run%d.txt",i))[,1]
}

# visualization
matplot(pde.t, pde.y, type = "l")
lines(exp.t,exp.y,lwd=2)

# remove last run due to numerical issue 
# with negative value
apply(pde.y, 2, function(x) sum(x<0))
rm.idx <- c(50)

# compute mse
# read table
mse <- read.table("../results/maxpro/maxpro_design_mse.txt")[,1]
#mse <- rep(0,nrow(D))
#for (i in 1:nrow(D)){
#  print(i)
#  if (i %in% rm.idx){
#    mse[i] <- Inf
#  } else {
#    exp.y.idw <- apply(matrix(exp.t,ncol=1), 1, idw, D=pde.t, y=pde.y[,i])
#    mse[i] <- mean((exp.y.idw-exp.y)^2)
#  }
#  print(mse[i])
#}
#write.table(mse, "../results/maxpro/maxpro_design_mse.txt", row.names = F, col.names = F)

# visualization
plot(exp.t, exp.y, "l")
lines(pde.t, pde.y[,which.min(mse)], col = "red", lty = 2)

# set random seed
set.seed(20210612)

# quadrature points for chi-square approximation
quad.node <- gauss.quad(10, kind = "legendre")
quad.node$nodes <- (quad.node$nodes+1)/2 # rescale to [0,1]
quad.node$weights <- quad.node$weights/2 # adjust the weight

# generalized chisq sp approximation
gchisq.sp <- NULL

# fpca setting
basis.breaks <- c(seq(0,5000,length.out=51)[1:50],seq(5000,120000,length.out=101))
basisobj <- create.bspline.basis(c(0,120000), breaks=basis.breaks)
fdParobj <- fdPar(basisobj, lambda = 2)

while (TRUE){
  print(nrow(D)+1)
  
  # functional pca
  print("Compute Functional PCA...")
  basis <- smooth.basis(pde.t, pde.y[,-rm.idx], fdParobj)
  fpca <- pca.fd(basis$fd, basisobj$nbasis)
  var.pca <- cumsum(fpca$varprop)
  nharm <- sum(var.pca < 0.99) + 1
  fpca <- pca.fd(basis$fd, nharm)
  # obtain eigenvalue
  fpca.eig.value <- fpca$values[1:nharm]
  # visualization
  # obtain mean and eigenvector
  # fpca.mean <- eval.fd(pde.t, fpca$meanfd)
  # fpca.eig.vector <- eval.fd(pde.t, fpca$harmonics)
  #for (idx in 1:nrow(D)){
  #  pde.pred <- fpca.mean
  #  for (i in 1:nharm){
  #    pde.pred <- pde.pred + fpca$scores[idx,i] * fpca.eig.vector[,i]
  #  }
  #  plot(pde.t, pde.y[,idx], "l", main = idx)
  #  lines(pde.t, pde.pred, lty = 2, col = "red")
  #}
  # obtain standardized pc score (mean 0, var 1)
  fpca.std.score <- t(t(fpca$scores)/sqrt(fpca.eig.value))
  # compute target standardized pc score
  fpca.mean <- eval.fd(exp.t, fpca$meanfd)
  target.basis <- smooth.basis(exp.t, matrix(exp.y-fpca.mean,ncol=1), fdParobj)
  target.score <- inprod(target.basis$fd, fpca$harmonics)
  target.std.score <- target.score / sqrt(fpca.eig.value)
  # visualization
  # fpca.eig.vector <- eval.fd(pde.t, fpca$harmonics)
  # target.pred <- fpca.mean
  # for (i in 1:nharm){
  # target.pred <- target.pred + target.score[i] * fpca.eig.vector[,i]
  #}
  #plot(exp.t, exp.y, "l")
  #lines(pde.t, target.pred, col = "red", lty = 2)
  
  # fit surrogate model
  print("Fit Surrogate Model...")
  models <- list()
  for (i in 1:nharm){
    invisible(capture.output(
      models[[i]] <- mlegp(D[-rm.idx,], fpca.std.score[,i])
    ))
  } 
  
  # obtain support points for the generalized chisq approximation
  if (is.null(gchisq.sp)){
    gchisq.sp <- sp(1000,nharm,dist.str=rep("normal",nharm))$sp
  } else{
    if (ncol(gchisq.sp)!=nharm){
      gchisq.sp <- sp(1000,nharm,dist.str=rep("normal",nharm))$sp
    }
  }
  
  # expected improvement for adaptive sampling
  print("Find Next Sample by EI...")
  mse.best <- min(mse)
  EI <- function(x){
    x <- exp(x) / (1+exp(x))
    score.fit <- rep(0, nharm)
    score.se <- rep(0, nharm)
    for (i in 1:nharm){
      score <- predict(models[[i]], x, se.fit = T)
      score.fit[i] <- score$fit
      score.se[i] <- score$se.fit
    }
    if (any(score.se < 1e-9)) return (0)
    # compute the expected improvement
    delta <- (score.fit - target.std.score)^2 / (score.se)^2
    lambda <- fpca.eig.value * (score.se)^2 / mse.best
    # integral approximation by quadrature
    quad.node.y <- rep(NA, length(quad.node$nodes))
    for (i in 1:length(quad.node$nodes)){
      quad.node.y[i] <- pgchisq(q=quad.node$nodes[i], lambda=lambda, delta=delta, D=gchisq.sp)
      # quad.node.y[i] <- 1 - imhof(q=quad.node$nodes[i], lambda=lambda, delta=delta)$Qq
    }
    val <- sum(quad.node.y * quad.node$weights) * mse.best
    return (-val) # negative for minimization
  }
  
  sobol.pts <- sobol(250, ncol(D), scrambling=1, seed=sample(1e9,1))
  sobol.pts[sobol.pts < 1e-3] <- 1e-3
  sobol.pts[sobol.pts > (1-1e-3)] <- (1-1e-3)
  sobol.pts <- log(sobol.pts/(1-sobol.pts))
  sobol.ei <- rep(0, nrow(sobol.pts))
  for (i in 1:nrow(sobol.pts)){
    sobol.ei[i] <- EI(sobol.pts[i,])
  }
  sobol.best <- sobol.pts[which.min(sobol.ei),] # negative ei
  x.cand <- optim(sobol.best, EI, method = "Nelder-Mead")
  x.next <- x.cand$par
  x.next <- exp(x.next) / (1 + exp(x.next))
  ei.best <- -x.cand$value # switch the sign due to minimization
  print(sprintf("ei: %.2e", ei.best))
  if (ei.best < 0){
    print("Cannot find a better point by EI!")
    break
  } else {
    print("Obtain the PDE output for the new sample point:")
    # transformation
    x.next.o <- x.next
    # Diffusivity, (1e-12,1e-9), log
    x.next.o[1] <- exp(log(1e-12) + x.next[1] * (log(1e-9) - log(1e-12)))
    # Surface Concentration, (4e-3,5e-3), original
    x.next.o[2] <- 4e-3 + x.next[2] * (5e-3 - 4e-3)
    # Accessible Polymer Concentration: (5e-3,6e-3), original
    x.next.o[3] <- 5e-3 + x.next[3] * (6e-3 - 5e-3)
    # Hindering, (500,2500), original
    x.next.o[4] <- 500 + x.next[4] * (2500 - 500)
    # Reaction rate, (1e-3,1e1), log
    x.next.o[5] <- exp(log(1e-3) + x.next[5] * (log(1e1) - log(1e-3)))
    print(x.next)
    print(x.next.o)
    write.table(c(nrow(D)+1, x.next.o), 
                file = "functional_new_design.txt", 
                row.names = F, col.names = F)
    # run matlab
    r <- tryCatch({
      system("/Applications/MATLAB_R2019a.app/bin/matlab -nodisplay -r followup", timeout = 900)
    }, warning = function(w){1})
    mass <- read.table(sprintf("../results/functional/run%d.txt",(nrow(D)+1)))[,1]
    pde.y <- cbind(pde.y, mass)
    # compute mse
    exp.y.idw <- apply(matrix(exp.t,ncol=1), 1, idw, D=pde.t, y=mass)
    mse <- c(mse, mean((exp.y.idw-exp.y)^2))
    write.table(mse, "../results/functional/all_design_mse.txt", row.names = F, col.names = F)
    # save all design
    D <- rbind(D, x.next)
    ini <- rbind(ini, c(nrow(D),x.next.o))
    write.csv(ini, "../results/functional/all_design.csv", row.names = F)
    # visualization
    plot(exp.t, exp.y, "l", main = nrow(D))
    # ground truth
    lines(pde.t, mass, col = "red", lty = 2)
    # gaussian process approximation
    score.fit <- rep(0, nharm)
    for (i in 1:nharm){
      score <- predict(models[[i]], x.next)
      score.fit[i] <- score
    }
    fpca.mean <- eval.fd(pde.t, fpca$meanfd)
    fpca.eig.vector <- eval.fd(pde.t, fpca$harmonics)
    pred <- fpca.mean
    for (i in 1:nharm){
      pred <- pred + score.fit[i] * sqrt(fpca.eig.value[i]) * fpca.eig.vector[,i]
    }
    lines(pde.t, pred, col = "blue", lty = 2)
  }
  
  print(sprintf("new mse: %.3f", mse[nrow(D)]))
  print(sprintf("min mse: %.3f", min(mse)))
  print("################################")
  
}

# delete functional new design files
unlink("functional_new_design.txt")

pairs(D, col=c(rep("red",50),rep("green",35)))

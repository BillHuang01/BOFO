rm(list=ls())

# load library
source("lib.R")
library(mlegp) # gaussian process
library(randtoolbox) # sobol point

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

# read computed mse
mse <- read.table("../results/maxpro/maxpro_design_mse.txt")[,1]

# visualization
plot(exp.t, exp.y, "l")
lines(pde.t, pde.y[,which.min(mse)], col = "red", lty = 2)

# set random seed
set.seed(20210612)

for (l in 1:30){
  print(nrow(D)+1)
  
  print("Fit Surrogate Model...")
  invisible(capture.output(model <- mlegp(D[-rm.idx,], mse[-rm.idx])))
  
  print("Find Next Sample by EI...")
  y.min <- min(mse)
  EI <- function(x){
    x <- exp(x) / (1 + exp(x)) # undo logit transform
    y.est <- predict(model, x, se.fit = T)
    y.fit <- y.est$fit
    y.se <- y.est$se.fit
    if (y.se < 1e-9) return (0)
    u <- (y.min - y.fit) / y.se
    val <- y.se * (u * pnorm(u) + dnorm(u))
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
                file = "standard_new_design.txt", 
                row.names = F, col.names = F)
    # run matlab
    r <- tryCatch({
      system("/Applications/MATLAB_R2019a.app/bin/matlab -nodisplay -r followup", timeout = 900)
    }, warning = function(w){1})
    mass <- read.table(sprintf("../results/standard/run%d.txt",(nrow(D)+1)))[,1]
    pde.y <- cbind(pde.y, mass)
    # compute mse
    exp.y.idw <- apply(matrix(exp.t,ncol=1), 1, idw, D=pde.t, y=mass)
    mse <- c(mse, mean((exp.y.idw-exp.y)^2))
    write.table(mse, "../results/standard/all_design_mse.txt", row.names = F, col.names = F)
    # save all design
    D <- rbind(D, x.next)
    ini <- rbind(ini, c(nrow(D),x.next.o))
    write.csv(ini, "../results/standard/all_design.csv", row.names = F)
    # visualization
    plot(exp.t, exp.y, "l", main = nrow(D))
    lines(pde.t, mass, col = "red", lty = 2)
  }
  
  print(sprintf("new mse: %.3f", mse[nrow(D)]))
  print(sprintf("min mse: %.3f", min(mse)))
  print("################################")
  
}

# delete standard new design files
unlink("standard_new_design.txt")

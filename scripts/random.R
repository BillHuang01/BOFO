rm(list=ls())

# load library
source("lib.R")

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
  
  # random sample a point from unit cube
  x.next <- runif(ncol(D))
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
              file = "random_new_design.txt", 
              row.names = F, col.names = F)
  # run matlab
  r <- tryCatch({
    system("/Applications/MATLAB_R2019a.app/bin/matlab -nodisplay -r followup", timeout = 900)
  }, warning = function(w){1})
  mass <- read.table(sprintf("../results/random/run%d.txt",(nrow(D)+1)))[,1]
  pde.y <- cbind(pde.y, mass)
  # compute mse
  exp.y.idw <- apply(matrix(exp.t,ncol=1), 1, idw, D=pde.t, y=mass)
  mse <- c(mse, mean((exp.y.idw-exp.y)^2))
  write.table(mse, "../results/random/all_design_mse.txt", row.names = F, col.names = F)
  # save all design
  D <- rbind(D, x.next)
  ini <- rbind(ini, c(nrow(D),x.next.o))
  write.csv(ini, "../results/random/all_design.csv", row.names = F)
  # visualization
  plot(exp.t, exp.y, "l", main = nrow(D))
  lines(pde.t, mass, col = "red", lty = 2)
  
  print(sprintf("new mse: %.3f", mse[nrow(D)]))
  print(sprintf("min mse: %.3f", min(mse)))
  print("################################")
  
}

# delete random new design files
unlink("random_new_design.txt")

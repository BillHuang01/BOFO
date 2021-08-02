# random seed
rm(list=ls())
set.seed(20210609)

# MaxPro Design Design 1
library(MaxPro)
p <- 5
D.mp <- MaxPro(MaxProLHD(p*10,p)$Design)$Design
pairs(D.mp)

# transform them to desired scale
D <- matrix(NA, nrow = p*10, ncol = p)
# Diffusivity, (1e-12,1e-9), log
D[,1] <- exp(log(1e-12) + D.mp[,1] * (log(1e-9) - log(1e-12)))
# Surface Concentration, (4e-3,5e-3), original
D[,2] <- 4e-3 + D.mp[,2] * (5e-3 - 4e-3)
# Accessible Polymer Concentration: (5e-3,6e-3), original
D[,3] <- 5e-3 + D.mp[,3] * (6e-3 - 5e-3)
# Hindering, (500,2500), original
D[,4] <- 500 + D.mp[,4] * (2500 - 500)
# Reaction rate, (1e-3,1e1), log
D[,5] <- exp(log(1e-3) + D.mp[,5] * (log(1e1) - log(1e-3)))
summary(D)

Design <- data.frame(
  run = c(1:nrow(D)),
  diffusivity = D[,1],
  surface_conc = D[,2],
  polymer_conc = D[,3],
  hindering = D[,4],
  reaction_rate = D[,5]
)
write.csv(Design, "data/maxpro_design.csv", row.names = F)

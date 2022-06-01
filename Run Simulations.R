# All Simulations
source("Functions.R")
library(knockoff)

# --------------------------------------------------
# Simulation Study 1
# --------------------------------------------------

# Knockoffs
sim <- read.csv('Simulation Study 1/Simulation_Study_Correlation_Knockoffs.csv')
names(sim)[1] <- 's'

args <- split(sim, seq(nrow(sim_var)))

for(row in 1:length(args)){
  do.call(apply_knockoffs, args[[row]])
}

# Lasso
sim <- read.csv('Simulation Study 1/Simulation_Study_Correlation_Lasso.csv')
names(sim)[1] <- 's'

args <- split(sim, seq(nrow(sim)))

for(row in 1:length(args)){
  do.call(apply_lasso, args[[row]])
}



# Adaptive Lasso
sim <- read.csv('Simulation Study 1/Simulation_Study_Correlation_Lasso_Adaptive.csv')
names(sim)[1] <- 's'

args <- split(sim, seq(nrow(sim)))

for(row in 1:length(args)){
  do.call(apply_lasso_adaptive, args[[row]])
}


# Elastic Net
sim <- read.csv('Simulation Study 1/Simulation_Study_Correlation_Elastic_Net.csv')
names(sim)[1] <- 's'

args <- split(sim, seq(nrow(sim)))

for(row in 1:length(args)){
  do.call(apply_elastic_net, args[[row]])
}

# ------------------------------------
# Simulation Study 2
# ------------------------------------

sim <- read.csv('Simulation Study 2/Simulation Study 2.csv')
names(sim)[1] <- 's'

args <- split(sim, seq(nrow(sim)))

for(row in 1:length(args)){
  do.call(apply_knockoffs_orp, args[[row]])
}


# --------------------------------------
# Simulation Study 3
# --------------------------------------
sim <- read.csv('Simulation Study 3/Simulation_Grid.csv')
names(sim)[1] <- 's'

args <- split(sim, seq(nrow(sim)))

for(row in 1:length(args)){
  do.call(apply_knockoffs, args[[row]])
}


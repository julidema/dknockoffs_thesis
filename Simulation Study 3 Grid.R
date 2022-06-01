# ---------------------------------
# Simulation Study 3 - Grid
# ---------------------------------
source("Functions.R")
library(knockoff)
library(ggplot2)
library(dplyr)

# Parameters of the runs
dfvals <- read.csv('Simulation Study 3/Simulation_Grid.csv')
# Keep n and p
dfvals <- dfvals[,c(2,3)]

# Results of the runs
res <- read.table('Simulation Study 3/fdpbin_grid.txt')
# Keep only fdp and power
res <- res[,c(1,2)]


# ---------------------------
# Data frame for plotting
# ---------------------------

df_plot <- data.frame(n = dfvals$n, p = dfvals$p, 
                      fdp = res$V1, power = res$V2)

# Compute mean fdp and power per n and p combination
df <- df_plot %>% group_by(n,p) %>% summarise(avfdp = mean(fdp), avpower = mean(power))



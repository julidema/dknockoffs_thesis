# ------------------------------------------------------------
# Simulation Study - Comparisons, Correlation and Predictions
# ------------------------------------------------------------
source("Functions.R")
library(ggplot2)
library(ggthemes)
library(wesanderson)

# -----------------------
# Knockoffs
# -----------------------

# FDR and Power
k_df_1 <- read.table('Simulation Study 1/Knockoffs/Binomial/fdpbin_corr00.txt')
k_df_2 <- read.table('Simulation Study 1/Knockoffs/Binomial/fdpbin_corr02.txt')
k_df_3 <- read.table('Simulation Study 1/Knockoffs/Binomial/fdpbin_corr04.txt')
k_df_4 <- read.table('Simulation Study 1/Knockoffs/Binomial/fdpbin_corr06.txt')
k_df_5 <- read.table('Simulation Study 1/Knockoffs/Binomial/fdpbin_corr08.txt')

k_fdr <- comp_mean_n(c(k_df_1$V1, k_df_2$V1, k_df_3$V1, k_df_4$V1, k_df_5$V1), 50)
k_power <- comp_mean_n(c(k_df_1$V2, k_df_2$V2, k_df_3$V2, k_df_4$V2, k_df_5$V2), 50)

# Number of Selected Variables
k_df_sel_1 <- read.table('Simulation Study 1/Knockoffs/Binomial/selbin_corr00.txt')
k_df_sel_2 <- read.table('Simulation Study 1/Knockoffs/Binomial/selbin_corr02.txt')
k_df_sel_3 <- read.table('Simulation Study 1/Knockoffs/Binomial/selbin_corr04.txt')
k_df_sel_4 <- read.table('Simulation Study 1/Knockoffs/Binomial/selbin_corr06.txt')
k_df_sel_5 <- read.table('Simulation Study 1/Knockoffs/Binomial/selbin_corr08.txt')

k_sel_1 <- comp_nr_selected(k_df_sel_1, 1)
k_sel_2 <- comp_nr_selected(k_df_sel_2, 1)
k_sel_3 <- comp_nr_selected(k_df_sel_3, 1)
k_sel_4 <- comp_nr_selected(k_df_sel_4, 1)
k_sel_5 <- comp_nr_selected(k_df_sel_5, 1)

k_sel_mean <- round(comp_mean_n(c(k_sel_1, k_sel_2, k_sel_3, k_sel_4, k_sel_5), 49))

# For density plots
d0_bin <- k_df_1[1:50, 2]
d1_bin <- k_df_2[1:50, 2]
d2_bin <- k_df_3[1:50, 2]
d3_bin <- k_df_4[1:50, 2]
d4_bin <- k_df_5[1:50, 2]
power_bin_dens <- c(d0_bin, d1_bin, d2_bin, d3_bin, d4_bin)


d0_bin <- k_df_1[1:50, 1]
d1_bin <- k_df_2[1:50, 1]
d2_bin <- k_df_3[1:50, 1]
d3_bin <- k_df_4[1:50, 1]
d4_bin <- k_df_5[1:50, 1]
fdp_bin_dens <- c(d0_bin, d1_bin, d2_bin, d3_bin, d4_bin)


# ---------------------------------
# Lasso
# ---------------------------------

l_df_1 <- read.table('Simulation Study 1/Lasso/Binomial/lasso_fdpbin_corr00.txt')
l_df_2 <- read.table('Simulation Study 1/Lasso/Binomial/lasso_fdpbin_corr02.txt')
l_df_3 <- read.table('Simulation Study 1/Lasso/Binomial/lasso_fdpbin_corr04.txt')
l_df_4 <- read.table('Simulation Study 1/Lasso/Binomial/lasso_fdpbin_corr06.txt')
l_df_5 <- read.table('Simulation Study 1/Lasso/Binomial/lasso_fdpbin_corr08.txt')

l_fdr <- comp_mean_n(c(l_df_1$V1, l_df_2$V1, l_df_3$V1, l_df_4$V1, l_df_5$V1), 50)
l_power <- comp_mean_n(c(l_df_1$V2, l_df_2$V2, l_df_3$V2, l_df_4$V2, l_df_5$V2), 50)


l_df_sel_1 <- read.table('Simulation Study 1/Lasso/Binomial/lasso_selbin_corr00.txt')
l_df_sel_2 <- read.table('Simulation Study 1/Lasso/Binomial/lasso_selbin_corr02.txt')
l_df_sel_3 <- read.table('Simulation Study 1/Lasso/Binomial/lasso_selbin_corr04.txt')
l_df_sel_4 <- read.table('Simulation Study 1/Lasso/Binomial/lasso_selbin_corr06.txt')
l_df_sel_5 <- read.table('Simulation Study 1/Lasso/Binomial/lasso_selbin_corr08.txt')

l_sel_1 <- comp_nr_selected(l_df_sel_1, 1)
l_sel_2 <- comp_nr_selected(l_df_sel_2, 1)
l_sel_3 <- comp_nr_selected(l_df_sel_3, 1)
l_sel_4 <- comp_nr_selected(l_df_sel_4, 1)
l_sel_5 <- comp_nr_selected(l_df_sel_5, 1)

l_sel_mean <- round(comp_mean_n(c(l_sel_1, l_sel_2, l_sel_3, l_sel_4, l_sel_5), 49))


# ----------------------------------
# Adaptive Lasso 
# ----------------------------------

al_df_1 <- read.table('Simulation Study 1/Adaptive Lasso/Binomial/adaptive_lasso_fdpbin_corr00.txt')
al_df_2 <- read.table('Simulation Study 1/Adaptive Lasso/Binomial/adaptive_lasso_fdpbin_corr02.txt')
al_df_3 <- read.table('Simulation Study 1/Adaptive Lasso/Binomial/adaptive_lasso_fdpbin_corr04.txt')
al_df_4 <- read.table('Simulation Study 1/Adaptive Lasso/Binomial/adaptive_lasso_fdpbin_corr06.txt')
al_df_5 <- read.table('Simulation Study 1/Adaptive Lasso/Binomial/adaptive_lasso_fdpbin_corr08.txt')

al_fdr <- comp_mean_n(c(al_df_1$V1, al_df_2$V1, al_df_3$V1, al_df_4$V1, al_df_5$V1), 50)
al_power <- comp_mean_n(c(al_df_1$V2, al_df_2$V2, al_df_3$V2, al_df_4$V2, al_df_5$V2), 50)


al_df_sel_1 <- read.table('Simulation Study 1/Adaptive Lasso/Binomial/adaptive_lasso_selbin_corr00.txt')
al_df_sel_2 <- read.table('Simulation Study 1/Adaptive Lasso/Binomial/adaptive_lasso_selbin_corr02.txt')
al_df_sel_3 <- read.table('Simulation Study 1/Adaptive Lasso/Binomial/adaptive_lasso_selbin_corr04.txt')
al_df_sel_4 <- read.table('Simulation Study 1/Adaptive Lasso/Binomial/adaptive_lasso_selbin_corr06.txt')
al_df_sel_5 <- read.table('Simulation Study 1/Adaptive Lasso/Binomial/adaptive_lasso_selbin_corr08.txt')

al_sel_1 <- comp_nr_selected(al_df_sel_1, 1)
al_sel_2 <- comp_nr_selected(al_df_sel_2, 1)
al_sel_3 <- comp_nr_selected(al_df_sel_3, 1)
al_sel_4 <- comp_nr_selected(al_df_sel_4, 1)
al_sel_5 <- comp_nr_selected(al_df_sel_5, 1)

al_sel_mean <- round(comp_mean_n(c(al_sel_1, al_sel_2, al_sel_3, al_sel_4, al_sel_5), 49))


# -------------------------
# Elastic Net
# -------------------------

en_df_1 <- read.table('Simulation Study 1/Elastic Net/Binomial/elastic_fdpbin_corr00.txt')
en_df_2 <- read.table('Simulation Study 1/Elastic Net/Binomial/elastic_fdpbin_corr02.txt')
en_df_3 <- read.table('Simulation Study 1/Elastic Net/Binomial/elastic_fdpbin_corr04.txt')
en_df_4 <- read.table('Simulation Study 1/Elastic Net/Binomial/elastic_fdpbin_corr06.txt')
en_df_5 <- read.table('Simulation Study 1/Elastic Net/Binomial/elastic_fdpbin_corr08.txt')

en_fdr <- comp_mean_n(c(en_df_1$V1, en_df_2$V1, en_df_3$V1, en_df_4$V1, en_df_5$V1), 50)
en_power <- comp_mean_n(c(en_df_1$V2, en_df_2$V2, en_df_3$V2, en_df_4$V2, en_df_5$V2), 50)


en_df_sel_1 <- read.table('Simulation Study 1/Elastic Net/Binomial/elastic_selbin_corr00.txt')
en_df_sel_2 <- read.table('Simulation Study 1/Elastic Net/Binomial/elastic_selbin_corr02.txt')
en_df_sel_3 <- read.table('Simulation Study 1/Elastic Net/Binomial/elastic_selbin_corr04.txt')
en_df_sel_4 <- read.table('Simulation Study 1/Elastic Net/Binomial/elastic_selbin_corr06.txt')
en_df_sel_5 <- read.table('Simulation Study 1/Elastic Net/Binomial/elastic_selbin_corr08.txt')

en_sel_1 <- comp_nr_selected(en_df_sel_1, 1)
en_sel_2 <- comp_nr_selected(en_df_sel_2, 1)
en_sel_3 <- comp_nr_selected(en_df_sel_3, 1)
en_sel_4 <- comp_nr_selected(en_df_sel_4, 1)
en_sel_5 <- comp_nr_selected(en_df_sel_5, 1)

en_sel_mean <- round(comp_mean_n(c(en_sel_1, en_sel_2, en_sel_3, en_sel_4, en_sel_5), 49))




# ----------------------------------
# Data frames for plotting
# ----------------------------------

sim_1_plotting <- data.frame(fdr = c(k_fdr, l_fdr, al_fdr, en_fdr),
                             power = c(k_power, l_power, al_power, en_power),
                             nsel = c(k_sel_mean, l_sel_mean, al_sel_mean, en_sel_mean),
                             Method = c(rep("Knockoffs+", 5), rep("LASSO", 5), rep("Adaptive LASSO", 5)
                                        , rep("Elastic Net", 5)),
                             corr = c(0, 0.2, 0.4, 0.6, 0.8))

sim_1_power <- data.frame(density = power_bin_dens, rho = factor(c(rep(0, 50), rep(0.2, 50), 
                                                                   rep(0.4, 50), rep(0.6,50),
                                                                   rep(0.8, 50))))


sim_1_fdp <- data.frame(density = fdp_bin_dens, rho = factor(c(rep(0, 50), rep(0.2, 50), 
                                                               rep(0.4, 50), rep(0.6,50),
                                                               rep(0.8, 50))))





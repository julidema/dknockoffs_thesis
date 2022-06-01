# -----------------------------------------------------------------------------------
# Simulation Study - Comparisons, Correlation and Predictions (Gaussin Response)
# -----------------------------------------------------------------------------------
source("Functions.R")
library(ggplot2)

# -----------------------
# Knockoffs
# -----------------------

# FDR and Power
k_df_1 <- read.table('Simulation Study 1/Knockoffs/Gaussian/fdpgaus_corr00.txt')
k_df_2 <- read.table('Simulation Study 1/Knockoffs/Gaussian/fdpgaus_corr02.txt')
k_df_3 <- read.table('Simulation Study 1/Knockoffs/Gaussian/fdpgaus_corr04.txt')
k_df_4 <- read.table('Simulation Study 1/Knockoffs/Gaussian/fdpgaus_corr06.txt')
k_df_5 <- read.table('Simulation Study 1/Knockoffs/Gaussian/fdpgaus_corr08.txt')

k_fdr <- comp_mean_n(c(k_df_1$V1, k_df_2$V1, k_df_3$V1, k_df_4$V1, k_df_5$V1), 50)
k_power <- comp_mean_n(c(k_df_1$V2, k_df_2$V2, k_df_3$V2, k_df_4$V2, k_df_5$V2), 50)

# Number of Selected Variables
k_df_sel_1 <- read.table('Simulation Study 1/Knockoffs/Gaussian/selgaus_corr00.txt')
k_df_sel_2 <- read.table('Simulation Study 1/Knockoffs/Gaussian/selgaus_corr02.txt')
k_df_sel_3 <- read.table('Simulation Study 1/Knockoffs/Gaussian/selgaus_corr04.txt')
k_df_sel_4 <- read.table('Simulation Study 1/Knockoffs/Gaussian/selgaus_corr06.txt')
k_df_sel_5 <- read.table('Simulation Study 1/Knockoffs/Gaussian/selgaus_corr08.txt')

k_sel_1 <- comp_nr_selected(k_df_sel_1, 1)
k_sel_2 <- comp_nr_selected(k_df_sel_2, 1)
k_sel_3 <- comp_nr_selected(k_df_sel_3, 1)
k_sel_4 <- comp_nr_selected(k_df_sel_4, 1)
k_sel_5 <- comp_nr_selected(k_df_sel_5, 1)


k_sel_mean <- round(comp_mean_n( c(k_sel_1, k_sel_2, k_sel_3, k_sel_4, k_sel_5), 49))

# For density plots
d0_gaus <- k_df_1[1:50, 2]
d1_gaus <- k_df_2[1:50, 2]
d2_gaus <- k_df_3[1:50, 2]
d3_gaus <- k_df_4[1:50, 2]
d4_gaus <- k_df_5[1:50, 2]
power_gaus_dens <- c(d0_gaus, d1_gaus, d2_gaus, d3_gaus, d4_gaus)


# ---------------------------------
# Lasso
# ---------------------------------

l_df_1 <- read.table('Simulation Study 1/Lasso/Gaussian/lasso_fdpgaus_corr00.txt')
l_df_2 <- read.table('Simulation Study 1/Lasso/Gaussian/lasso_fdpgaus_corr02.txt')
l_df_3 <- read.table('Simulation Study 1/Lasso/Gaussian/lasso_fdpgaus_corr04.txt')
l_df_4 <- read.table('Simulation Study 1/Lasso/Gaussian/lasso_fdpgaus_corr06.txt')
l_df_5 <- read.table('Simulation Study 1/Lasso/Gaussian/lasso_fdpgaus_corr08.txt')

l_fdr <- comp_mean_n(c(l_df_1$V1, l_df_2$V1, l_df_3$V1, l_df_4$V1, l_df_5$V1), 50)
l_power <- comp_mean_n(c(l_df_1$V2, l_df_2$V2, l_df_3$V2, l_df_4$V2, l_df_5$V2), 50)


l_df_sel_1 <- read.table('Simulation Study 1/Lasso/Gaussian/lasso_selgaus_corr00.txt')
l_df_sel_2 <- read.table('Simulation Study 1/Lasso/Gaussian/lasso_selgaus_corr02.txt')
l_df_sel_3 <- read.table('Simulation Study 1/Lasso/Gaussian/lasso_selgaus_corr04.txt')
l_df_sel_4 <- read.table('Simulation Study 1/Lasso/Gaussian/lasso_selgaus_corr06.txt')
l_df_sel_5 <- read.table('Simulation Study 1/Lasso/Gaussian/lasso_selgaus_corr08.txt')

l_sel_1 <- comp_nr_selected(l_df_sel_1, 1)
l_sel_2 <- comp_nr_selected(l_df_sel_2, 1)
l_sel_3 <- comp_nr_selected(l_df_sel_3, 1)
l_sel_4 <- comp_nr_selected(l_df_sel_4, 1)
l_sel_5 <- comp_nr_selected(l_df_sel_5, 1)


l_sel_mean <- round(comp_mean_n(c(l_sel_1, l_sel_2, l_sel_3, l_sel_4, l_sel_5), 49))


# ----------------------------------
# Adaptive Lasso 
# ----------------------------------

al_df_1 <- read.table('Simulation Study 1/Adaptive Lasso/Gaussian/adaptive_lasso_fdpgaus_corr00.txt')
al_df_2 <- read.table('Simulation Study 1/Adaptive Lasso/Gaussian/adaptive_lasso_fdpgaus_corr02.txt')
al_df_3 <- read.table('Simulation Study 1/Adaptive Lasso/Gaussian/adaptive_lasso_fdpgaus_corr04.txt')
al_df_4 <- read.table('Simulation Study 1/Adaptive Lasso/Gaussian/adaptive_lasso_fdpgaus_corr06.txt')
al_df_5 <- read.table('Simulation Study 1/Adaptive Lasso/Gaussian/adaptive_lasso_fdpgaus_corr08.txt')

al_fdr <- comp_mean_n(c(al_df_1$V1, al_df_2$V1, al_df_3$V1, al_df_4$V1, al_df_5$V1), 50)
al_power <- comp_mean_n(c(al_df_1$V2, al_df_2$V2, al_df_3$V2, al_df_4$V2, al_df_5$V2), 50)


al_df_sel_1 <- read.table('Simulation Study 1/Adaptive Lasso/Gaussian/adaptive_lasso_selgaus_corr00.txt')
al_df_sel_2 <- read.table('Simulation Study 1/Adaptive Lasso/Gaussian/adaptive_lasso_selgaus_corr02.txt')
al_df_sel_3 <- read.table('Simulation Study 1/Adaptive Lasso/Gaussian/adaptive_lasso_selgaus_corr04.txt')
al_df_sel_4 <- read.table('Simulation Study 1/Adaptive Lasso/Gaussian/adaptive_lasso_selgaus_corr06.txt')
al_df_sel_5 <- read.table('Simulation Study 1/Adaptive Lasso/Gaussian/adaptive_lasso_selgaus_corr08.txt')

al_sel_1 <- comp_nr_selected(al_df_sel_1, 1)
al_sel_2 <- comp_nr_selected(al_df_sel_2, 1)
al_sel_3 <- comp_nr_selected(al_df_sel_3, 1)
al_sel_4 <- comp_nr_selected(al_df_sel_4, 1)
al_sel_5 <- comp_nr_selected(al_df_sel_5, 1)

al_sel_mean <- round(comp_mean_n(c(al_sel_1, al_sel_2, al_sel_3, al_sel_4, al_sel_5), 49))


# -------------------------
# Elastic Net
# -------------------------

en_df_1 <- read.table('Simulation Study 1/Elastic Net/Gaussian/elastic_fdpgaus_corr00.txt')
en_df_2 <- read.table('Simulation Study 1/Elastic Net/Gaussian/elastic_fdpgaus_corr02.txt')
en_df_3 <- read.table('Simulation Study 1/Elastic Net/Gaussian/elastic_fdpgaus_corr04.txt')
en_df_4 <- read.table('Simulation Study 1/Elastic Net/Gaussian/elastic_fdpgaus_corr06.txt')
en_df_5 <- read.table('Simulation Study 1/Elastic Net/Gaussian/elastic_fdpgaus_corr08.txt')

en_fdr <- comp_mean_n(c(en_df_1$V1, en_df_2$V1, en_df_3$V1, en_df_4$V1, en_df_5$V1), 50)
en_power <- comp_mean_n(c(en_df_1$V2, en_df_2$V2, en_df_3$V2, en_df_4$V2, en_df_5$V2), 50)


en_df_sel_1 <- read.table('Simulation Study 1/Elastic Net/Gaussian/elastic_selgaus_corr00.txt')
en_df_sel_2 <- read.table('Simulation Study 1/Elastic Net/Gaussian/elastic_selgaus_corr02.txt')
en_df_sel_3 <- read.table('Simulation Study 1/Elastic Net/Gaussian/elastic_selgaus_corr04.txt')
en_df_sel_4 <- read.table('Simulation Study 1/Elastic Net/Gaussian/elastic_selgaus_corr06.txt')
en_df_sel_5 <- read.table('Simulation Study 1/Elastic Net/Gaussian/elastic_selgaus_corr08.txt')

en_sel_1 <- comp_nr_selected(en_df_sel_1, 1)
en_sel_2 <- comp_nr_selected(en_df_sel_2, 1)
en_sel_3 <- comp_nr_selected(en_df_sel_3, 1)
en_sel_4 <- comp_nr_selected(en_df_sel_4, 1)
en_sel_5 <- comp_nr_selected(en_df_sel_5, 1)

en_sel_mean <- round(comp_mean_n(c(en_sel_1, en_sel_2, en_sel_3, en_sel_4, en_sel_5), 49))




# -----------------------------------------------
# Data Frames for plotting
# -----------------------------------------------

sim_1_gaus_plotting <- data.frame(fdr = c(k_fdr, l_fdr, al_fdr, en_fdr),
                                  power = c(k_power, l_power, al_power, en_power),
                                  nsel = c(k_sel_mean, l_sel_mean, al_sel_mean, en_sel_mean),
                                  Method = c(rep("Knockoffs+", 5), rep("LASSO", 5), rep("Adaptive LASSO", 5),
                                             rep("Elastic Net", 5)),
                                  corr = c(0, 0.2, 0.4, 0.6, 0.8))


sim_1_power_gaus <- data.frame(density = power_gaus_dens, rho = factor(c(rep(0, 50), rep(0.2, 50), 
                                                                   rep(0.4, 50), rep(0.6,50),
                                                                   rep(0.8, 50))))

sum(k_sel_1 == 0)
sum(k_sel_2 == 0)
sum(k_sel_3 == 0)
sum(k_sel_4 == 0)
sum(k_sel_5 == 0)



# Lasso
sum(l_sel_1 == 0)
sum(l_sel_2 == 0)
sum(l_sel_3 == 0)
sum(l_sel_4 == 0)
sum(l_sel_5 == 0)


# Adaptive Lasso
sum(al_sel_1 == 0)
sum(al_sel_2 == 0)
sum(al_sel_3 == 0)
sum(al_sel_4 == 0)
sum(al_sel_5 == 0)


# Elastic Net
sum(en_sel_1 == 0)
sum(en_sel_2 == 0)
sum(en_sel_3 == 0)
sum(en_sel_4 == 0)
sum(en_sel_5 == 0)



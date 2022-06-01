# -----------------------------------------------------------------------------------
# Simulation Study - Comparisons, Correlation and Predictions (Poisson Response)
# -----------------------------------------------------------------------------------
source("Functions.R")
library(ggplot2)
library(dplyr)

# -----------------------
# Knockoffs
# -----------------------

# FDR and Power
k_df_1 <- read.table('Simulation Study 1/Knockoffs/Poisson/fdppois_corr00.txt')
k_df_2 <- read.table('Simulation Study 1/Knockoffs/Poisson/fdppois_corr02.txt')
k_df_3 <- read.table('Simulation Study 1/Knockoffs/Poisson/fdppois_corr04.txt')
k_df_4 <- read.table('Simulation Study 1/Knockoffs/Poisson/fdppois_corr06.txt')
k_df_5 <- read.table('Simulation Study 1/Knockoffs/Poisson/fdppois_corr08.txt')

k_fdr <- comp_mean_n(c(k_df_1$V1, k_df_2$V1, k_df_3$V1, k_df_4$V1, k_df_5$V1), 50)
k_power <- comp_mean_n(c(k_df_1$V2, k_df_2$V2, k_df_3$V2, k_df_4$V2, k_df_5$V2), 50)

# Number of Selected Variables
k_df_sel_1 <- read.table('Simulation Study 1/Knockoffs/Poisson/selpois_corr00.txt')
k_df_sel_2 <- read.table('Simulation Study 1/Knockoffs/Poisson/selpois_corr02.txt')
k_df_sel_3 <- read.table('Simulation Study 1/Knockoffs/Poisson/selpois_corr04.txt')
k_df_sel_4 <- read.table('Simulation Study 1/Knockoffs/Poisson/selpois_corr06.txt')
k_df_sel_5 <- read.table('Simulation Study 1/Knockoffs/Poisson/selpois_corr08.txt')

k_sel_1 <- comp_nr_selected(k_df_sel_1, 1)
k_sel_2 <- comp_nr_selected(k_df_sel_2, 1)
k_sel_3 <- comp_nr_selected(k_df_sel_3, 1)
k_sel_4 <- comp_nr_selected(k_df_sel_4, 1)
k_sel_5 <- comp_nr_selected(k_df_sel_5, 1)


k_sel_mean <- round(comp_mean_n( c(k_sel_1, k_sel_2, k_sel_3, k_sel_4, k_sel_5), 49))


d0_pois <- k_df_1[1:50, 2]
d1_pois <- k_df_2[1:50, 2]
d2_pois <- k_df_3[1:50, 2]
d3_pois <- k_df_4[1:50, 2]
d4_pois <- k_df_5[1:50, 2]
power_pois_dens <- c(d0_pois, d1_pois, d2_pois, d3_pois, d4_pois)




# ---------------------------------
# Lasso
# ---------------------------------

l_df_1 <- read.table('Simulation Study 1/Lasso/Poisson/lasso_fdppois_corr00.txt')
l_df_2 <- read.table('Simulation Study 1/Lasso/Poisson/lasso_fdppois_corr02.txt')
l_df_3 <- read.table('Simulation Study 1/Lasso/Poisson/lasso_fdppois_corr04.txt')
l_df_4 <- read.table('Simulation Study 1/Lasso/Poisson/lasso_fdppois_corr06.txt')
l_df_5 <- read.table('Simulation Study 1/Lasso/Poisson/lasso_fdppois_corr08.txt')

l_fdr <- comp_mean_n(c(l_df_1$V1, l_df_2$V1, l_df_3$V1, l_df_4$V1, l_df_5$V1), 50)
l_power <- comp_mean_n(c(l_df_1$V2, l_df_2$V2, l_df_3$V2, l_df_4$V2, l_df_5$V2), 50)


l_df_sel_1 <- read.table('Simulation Study 1/Lasso/Poisson/lasso_selpois_corr00.txt')
l_df_sel_2 <- read.table('Simulation Study 1/Lasso/Poisson/lasso_selpois_corr02.txt')
l_df_sel_3 <- read.table('Simulation Study 1/Lasso/Poisson/lasso_selpois_corr04.txt')
l_df_sel_4 <- read.table('Simulation Study 1/Lasso/Poisson/lasso_selpois_corr06.txt')
l_df_sel_5 <- read.table('Simulation Study 1/Lasso/Poisson/lasso_selpois_corr08.txt')

l_sel_1 <- comp_nr_selected(l_df_sel_1, 1)
l_sel_2 <- comp_nr_selected(l_df_sel_2, 1)
l_sel_3 <- comp_nr_selected(l_df_sel_3, 1)
l_sel_4 <- comp_nr_selected(l_df_sel_4, 1)
l_sel_5 <- comp_nr_selected(l_df_sel_5, 1)

l_sel_1_rm <- mean(l_sel_1[l_sel_1 != 0])
l_sel_2_rm <- mean(l_sel_2[l_sel_2 != 0])
l_sel_3_rm <- mean(l_sel_3[l_sel_3 != 0])
l_sel_4_rm <- mean(l_sel_4[l_sel_4 != 0])
l_sel_5_rm <- mean(l_sel_5[l_sel_5 != 0])


l_sel_mean <- round(comp_mean_n(c(l_sel_1, l_sel_2, l_sel_3, l_sel_4, l_sel_5), 49))
l_sel_mean_rm <- round(c(l_sel_1_rm, l_sel_2_rm, l_sel_3_rm, l_sel_4_rm, l_sel_5_rm))


l_fp_rm <- data.frame(fdr = c(l_df_1$V1, l_df_2$V1, l_df_3$V1, l_df_4$V1, l_df_5$V1),
                      power = c(l_df_1$V2, l_df_2$V2, l_df_3$V2, l_df_4$V2, l_df_5$V2),
                      sel = c(l_sel_1, l_sel_2, l_sel_3, l_sel_4, l_sel_5),
                      corr = c(rep(0, 50), rep(0.2, 50), rep(0.4, 50), rep(0.6, 50), rep(0.8, 50)),
                      Method = "LASSO")

# ----------------------------------
# Adaptive Lasso
# ----------------------------------

al_df_1 <- read.table('Simulation Study 1/Adaptive Lasso/Poisson/adaptive_lasso_fdppois_corr00.txt')
al_df_2 <- read.table('Simulation Study 1/Adaptive Lasso/Poisson/adaptive_lasso_fdppois_corr02.txt')
al_df_3 <- read.table('Simulation Study 1/Adaptive Lasso/Poisson/adaptive_lasso_fdppois_corr04.txt')
al_df_4 <- read.table('Simulation Study 1/Adaptive Lasso/Poisson/adaptive_lasso_fdppois_corr06.txt')
al_df_5 <- read.table('Simulation Study 1/Adaptive Lasso/Poisson/adaptive_lasso_fdppois_corr08.txt')

al_fdr <- comp_mean_n(c(al_df_1$V1, al_df_2$V1, al_df_3$V1, al_df_4$V1, al_df_5$V1), 50)
al_power <- comp_mean_n(c(al_df_1$V2, al_df_2$V2, al_df_3$V2, al_df_4$V2, al_df_5$V2), 50)


al_df_sel_1 <- read.table('Simulation Study 1/Adaptive Lasso/Poisson/adaptive_lasso_selpois_corr00.txt')
al_df_sel_2 <- read.table('Simulation Study 1/Adaptive Lasso/Poisson/adaptive_lasso_selpois_corr02.txt')
al_df_sel_3 <- read.table('Simulation Study 1/Adaptive Lasso/Poisson/adaptive_lasso_selpois_corr04.txt')
al_df_sel_4 <- read.table('Simulation Study 1/Adaptive Lasso/Poisson/adaptive_lasso_selpois_corr06.txt')
al_df_sel_5 <- read.table('Simulation Study 1/Adaptive Lasso/Poisson/adaptive_lasso_selpois_corr08.txt')

al_sel_1 <- comp_nr_selected(al_df_sel_1, 1)
al_sel_2 <- comp_nr_selected(al_df_sel_2, 1)
al_sel_3 <- comp_nr_selected(al_df_sel_3, 1)
al_sel_4 <- comp_nr_selected(al_df_sel_4, 1)
al_sel_5 <- comp_nr_selected(al_df_sel_5, 1)

al_sel_mean <- round(comp_mean_n(c(al_sel_1, al_sel_2, al_sel_3, al_sel_4, al_sel_5), 49))

al_sel_1_rm <- mean(al_sel_1[al_sel_1 != 0])
al_sel_2_rm <- mean(al_sel_2[al_sel_2 != 0])
al_sel_3_rm <- mean(al_sel_3[al_sel_3 != 0])
al_sel_4_rm <- mean(al_sel_4[al_sel_4 != 0])
al_sel_5_rm <- mean(al_sel_5[al_sel_5 != 0])
al_sel_mean_rm <- round(c(al_sel_1_rm, al_sel_2_rm, al_sel_3_rm, al_sel_4_rm, al_sel_5_rm))


al_fp_rm <- data.frame(fdr = c(al_df_1$V1, al_df_2$V1, al_df_3$V1, al_df_4$V1, al_df_5$V1),
                       power = c(al_df_1$V2, al_df_2$V2, al_df_3$V2, al_df_4$V2, al_df_5$V2),
                       sel = c(al_sel_1, al_sel_2, al_sel_3, al_sel_4, al_sel_5),
                       corr = c(rep(0, 50), rep(0.2, 50), rep(0.4, 50), rep(0.6, 50), rep(0.8, 50)),
                       Method = "Adaptive LASSO")
# -------------------------
# Elastic Net
# -------------------------

en_df_1 <- read.table('Simulation Study 1/Elastic Net/Poisson/elastic_fdppois_corr00.txt')
en_df_2 <- read.table('Simulation Study 1/Elastic Net/Poisson/elastic_fdppois_corr02.txt')
en_df_3 <- read.table('Simulation Study 1/Elastic Net/Poisson/elastic_fdppois_corr04.txt')
en_df_4 <- read.table('Simulation Study 1/Elastic Net/Poisson/elastic_fdppois_corr06.txt')
en_df_5 <- read.table('Simulation Study 1/Elastic Net/Poisson/elastic_fdppois_corr08.txt')

en_fdr <- comp_mean_n(c(en_df_1$V1, en_df_2$V1, en_df_3$V1, en_df_4$V1, en_df_5$V1), 50)
en_power <- comp_mean_n(c(en_df_1$V2, en_df_2$V2, en_df_3$V2, en_df_4$V2, en_df_5$V2), 50)


en_df_sel_1 <- read.table('Simulation Study 1/Elastic Net/Poisson/elastic_selpois_corr00.txt')
en_df_sel_2 <- read.table('Simulation Study 1/Elastic Net/Poisson/elastic_selpois_corr02.txt')
en_df_sel_3 <- read.table('Simulation Study 1/Elastic Net/Poisson/elastic_selpois_corr04.txt')
en_df_sel_4 <- read.table('Simulation Study 1/Elastic Net/Poisson/elastic_selpois_corr06.txt')
en_df_sel_5 <- read.table('Simulation Study 1/Elastic Net/Poisson/elastic_selpois_corr08.txt')

# Average number of selected variables
en_sel_1 <- comp_nr_selected(en_df_sel_1, 1)
en_sel_2 <- comp_nr_selected(en_df_sel_2, 1)
en_sel_3 <- comp_nr_selected(en_df_sel_3, 1)
en_sel_4 <- comp_nr_selected(en_df_sel_4, 1)
en_sel_5 <- comp_nr_selected(en_df_sel_5, 1)

en_sel_mean <- round(comp_mean_n(c(en_sel_1, en_sel_2, en_sel_3, en_sel_4, en_sel_5), 49))


# Average number of selected variables removing no selection cases
en_sel_1_rm <- mean(en_sel_1[en_sel_1 != 0])
en_sel_2_rm <- mean(en_sel_2[en_sel_2 != 0])
en_sel_3_rm <- mean(en_sel_3[en_sel_3 != 0])
en_sel_4_rm <- mean(en_sel_4[en_sel_4 != 0])
en_sel_5_rm <- mean(en_sel_5[en_sel_5 != 0])
en_sel_mean_rm <- round(c(en_sel_1_rm, en_sel_2_rm, en_sel_3_rm, en_sel_4_rm, en_sel_5_rm))


en_fp_rm <- data.frame(fdr = c(en_df_1$V1, en_df_2$V1, en_df_3$V1, en_df_4$V1, en_df_5$V1),
                       power = c(en_df_1$V2, en_df_2$V2, en_df_3$V2, en_df_4$V2, en_df_5$V2),
                       sel = c(en_sel_1, en_sel_2, en_sel_3, en_sel_4, en_sel_5),
                       corr = c(rep(0, 50), rep(0.2, 50), rep(0.4, 50), rep(0.6, 50), rep(0.8, 50)),
                       Method = "Elastic Net")


# -------------------------------
# Data Frames for plotting
# -------------------------------

sim_1_pois_plotting <- data.frame(fdr = c(k_fdr, l_fdr, al_fdr, en_fdr),
                                  power = c(k_power, l_power, al_power, en_power),
                                  nsel = c(k_sel_mean, l_sel_mean, al_sel_mean, en_sel_mean),
                                  nselrm = c(k_sel_mean, l_sel_mean_rm, al_sel_mean_rm, en_sel_mean_rm),
                                  Method = c(rep("Knockoffs+", 5), rep("LASSO", 5), rep("Adaptive LASSO", 5),
                                             rep("Elastic Net", 5)),
                                  corr = c(0, 0.2, 0.4, 0.6, 0.8))



# FOR TABLE 1
# Look at the number of times a selection is performed vs no selection perormed
# And then look at the average number of selected variables when a selection is performed

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


# For plotting fdr and power with zero selection entries removed
fdp_power_nozero <- rbind(l_fp_rm, al_fp_rm, en_fp_rm)

nozero <- fdp_power_nozero %>% filter(sel != 0) %>%
  group_by(Method, corr) %>% summarise(avfdp = mean(fdr),
                                       avpower = mean(power))




# ---------------------------------------------------------------------------
# Simulation Study Comparing Knockoff to Knockoff+ while varying amplitude
# ---------------------------------------------------------------------------
source("Functions.R")
library(ggplot2)

# ------------
# Knockoffs
# ------------

# FDP and Power
# A = 1
df_1 <- read.table("Simulation Study 2/fdpbin_a1_q1.txt")
# A = 5
df_5 <- read.table("Simulation Study 2/fdpbin_a5_q1.txt")
# A = 10
df_10 <- read.table("Simulation Study 2/fdpbin_a10_q1.txt")
# A = 20
df_20 <- read.table("Simulation Study 2/fdpbin_a20_q1.txt")
# A = 30
df_30 <- read.table("Simulation Study 2/fdpbin_a30_q1.txt")
# A = 40
df_40 <- read.table("Simulation Study 2/fdpbin_a40_q1.txt")
# A = 50
df_50 <- read.table("Simulation Study 2/fdpbin_a50_q1.txt")


# FDP and Power density
df_1 <- df_1[,c("V1", "V2")]
df_5 <- df_5[,c("V1", "V2")]
df_10 <- df_10[,c("V1", "V2")]
df_20 <- df_20[,c("V1", "V2")]
df_30 <- df_30[,c("V1", "V2")]
df_40 <- df_40[,c("V1", "V2")]
df_50 <- df_50[,c("V1", "V2")]



# --------------------------
# Knockoffs+
# --------------------------

# FDP and Power
# A = 1
df_p1 <- read.table("Simulation Study 2/fdpbin_pa1_q1.txt")
# A = 5
df_p5 <- read.table("Simulation Study 2/fdpbin_pa5_q1.txt")
# A = 10
df_p10 <- read.table("Simulation Study 2/fdpbin_pa10_q1.txt")
# A = 20
df_p20 <- read.table("Simulation Study 2/fdpbin_pa20_q1.txt")
# A = 30
df_p30 <- read.table("Simulation Study 2/fdpbin_pa30_q1.txt")
# A = 40
df_p40 <- read.table("Simulation Study 2/fdpbin_pa40_q1.txt")
# A = 50
df_p50 <- read.table("Simulation Study 2/fdpbin_pa50_q1.txt")


# FDP and Power density
df_p1 <- df_p1[,c("V1", "V2")]
df_p5 <- df_p5[,c("V1", "V2")]
df_p10 <- df_p10[,c("V1", "V2")]
df_p20 <- df_p20[,c("V1", "V2")]
df_p30 <- df_p30[,c("V1", "V2")]
df_p40 <- df_p40[,c("V1", "V2")]
df_p50 <- df_p50[,c("V1", "V2")]


# -----------------------
# Both methods
# -----------------------


# ALL FDP and Power

knockoffs <- rbind(df_1, df_5, df_10, df_20, df_30, df_40, df_50)
knockoffs_p <- rbind(df_p1, df_p5, df_p10, df_p20, df_p30, df_p40, df_p50)


# Average FDP and power per amplitude 
fdr <- comp_mean_n(knockoffs$V1, 20)
fdrp <- comp_mean_n(knockoffs_p$V1, 20)
power <- comp_mean_n(knockoffs$V2, 20)
powerp <- comp_mean_n(knockoffs_p$V2, 20)

# ---------------------------
# Data Frames for plotting
# ---------------------------

# Average FDP and Power
sim_plus <- data.frame(fdr = c(fdr, fdrp), power = c(power, powerp), 
                       amp = c(rep(1,5), rep(5,5), rep(10,5), rep(20,5), rep(30,5), rep(40,5), rep(50,5)),
                       Method = c(rep("Knockoffs", 35), rep("Knockoffs+", 35)),
                       q = as.factor(c(0.05, 0.1, 0.2, 0.3, 0.5)),
                       qline = c(0.05, 0.1, 0.2, 0.3, 0.5))








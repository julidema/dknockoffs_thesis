# ------------------------------------------------------------
# Simulation Study - Effect of Ratio of p/n on Knockoff Procedure 
# ------------------------------------------------------------
source("Functions.R")
library(ggplot2)


# ------------------------------
# n Fixed
# ------------------------------
ratio_n <- c(0.01, seq(0.05, 0.5, by = 0.05), seq(0.6, 1, by = 0.1))

# Parameters of the simulaitons
df_n_500 <- read.csv("Ratio FINAL/Simulation_Study_n500_FINAL.csv")[1:528,]
df_n_1000 <- read.csv("Ratio FINAL/Simulation_Study_n1000_FINAL.csv")[1:528,]
df_n_2000 <- read.csv("Ratio FINAL/Simulation_Study_n2000_FINAL.csv")[1:528,]

df_n_500 <- df_n_500[,c('n','p','k')]
df_n_1000 <- df_n_1000[,c('n','p','k')]
df_n_2000 <- df_n_2000[,c('n','p','k')]


# Results of the simulations
res_n_500 <- read.table('Ratio FINAL/fdpbin_nfix500.txt')
res_n_1000 <- read.table('Ratio FINAL/fdpbin_nfix1000.txt')
res_n_2000 <- read.table('Ratio FINAL/fdpbin_nfix2000.txt')

res_n_500 <- res_n_500[, c('V1', 'V2')][1:528,]
res_n_1000 <- res_n_1000[, c('V1', 'V2')][1:528,]
res_n_2000 <- res_n_2000[, c('V1', 'V2')][1:528,]


n_fix <- rbind(cbind(df_n_500, res_n_500), cbind(df_n_1000, res_n_1000), 
               cbind(df_n_2000, res_n_2000))

n_fix[is.na(n_fix)] <- 0
dim(n_fix)

test <- rep(rep(ratio_n, each = 50), 2)[1:1584]

n_fix$ratio <- rep(ratio_n, each = 50)
n_fix$ratio <- test
names(n_fix)[4:5] <- c("fdp", 'power')


n_fix_dropp <- n_fix[,-2]

n_fix_sum <- n_fix_dropp %>% group_by(n, ratio) %>% summarise(avfdp = mean(fdp), 
                                                            avpower = mean(power))

n_fixdf <- as.data.frame(n_fix_sum)
n_fixdf$n <- as.factor(n_fixdf$n)

ggplot(n_fixdf, aes(x = ratio)) +
  geom_line(aes(y = avpower, color = n), size = 1.5) + 
  geom_hline(yintercept = 0.1, linetype = 2) + 
  labs(y = 'Power', x = "Ratio p/n") +
  theme_minimal() + 
  guides(color = guide_legend(title = "Number of observations n (Fixed)"))+
  theme(legend.position = c(0.8,0.8), 
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_text("Number of covariates", face = "bold"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))+ 
  # scale_x_reverse()+
  scale_color_brewer(palette = "Set2")




ggplot(pfixdf, aes(x = ratio)) +
  geom_line(aes(y = avfdp, color = p), size = 1.5) + 
  geom_hline(yintercept = 0.1, linetype = 2) + 
  labs(y = 'Power', x = "Ratio p/n") +
  theme_minimal() + 
  guides(color = guide_legend(title = "Number of covariates p (Fixed)"))+
  theme(legend.position = c(0.8,0.8), 
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_text("Number of covariates", face = "bold"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))+ 
  # scale_x_reverse()+
  scale_color_brewer(palette = "Set2")


# ----------------------------
# p Fixed
# ----------------------------

ratio_p <- c(0.01, seq(0.05, 0.5, by = 0.05), seq(0.6, 1, by = 0.1), 1.5, 2)

# Parameters of the simulations
df_p_100 <- read.csv('Ratio FINAL/Simulation_Study_p100_FINAL.csv')
df_p_200 <- read.csv('Ratio FINAL/Simulation_Study_p200_FINAL.csv')
df_p_300 <- read.csv('Ratio FINAL/Simulation_Study_p300_FINAL.csv')

df_p_100 <- df_p_100[,c('n','p','k')]
df_p_200 <- df_p_200[,c('n','p','k')]
df_p_300 <- df_p_300[,c('n','p','k')]


# Results of the simulations
res_p_100 <- read.table("Ratio FINAL/fdpbin_pfix100.txt")
res_p_200 <- read.table("Ratio FINAL/fdpbin_pfix200.txt")
res_p_300 <- read.table("Ratio FINAL/fdpbin_pfix300.txt")


res_p_100 <- res_p_100[, c('V1', 'V2')]
res_p_200 <- res_p_200[, c('V1', 'V2')]
res_p_300 <- res_p_300[, c('V1', 'V2')]


p_fix <- rbind(cbind(df_p_100, res_p_100), cbind(df_p_200, res_p_200), cbind(df_p_300, res_p_300))
p_fix$ratio <- rep(ratio_p, each = 50)
names(p_fix)[4:5] <- c("fdp", 'power')


test <- comp_mean_n(p_fix$power, 50)
length(test)



pfix_sum <- pfix_dropn %>% group_by(p, ratio) %>% summarise(avfdp = mean(fdp), 
                                                   avpower = mean(power))
head(pfix_sum)

pfixdf <- as.data.frame(pfix_sum)
pfixdf$p <- as.factor(pfixdf$p)

head(pfixdf)

ggplot(pfixdf, aes(x = ratio)) + 
  geom_line(aes(y = avpower, color = p))


ggplot(pfixdf, aes(x = ratio)) +
  geom_line(aes(y = avpower, color = p), size = 1.5) + 
  geom_hline(yintercept = 0.1, linetype = 2) + 
  labs(y = 'Power', x = "Ratio p/n") +
  theme_minimal() + 
  guides(color = guide_legend(title = "Number of covariates p (Fixed)"))+
  theme(legend.position = c(0.8,0.8), 
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_text("Number of covariates", face = "bold"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))+ 
  # scale_x_reverse()+
  scale_color_brewer(palette = "Set2")




ggplot(pfixdf, aes(x = ratio)) +
  geom_line(aes(y = avfdp, color = p), size = 1.5) + 
  geom_hline(yintercept = 0.1, linetype = 2) + 
  labs(y = 'Power', x = "Ratio p/n") +
  theme_minimal() + 
  guides(color = guide_legend(title = "Number of covariates p (Fixed)"))+
  theme(legend.position = c(0.8,0.8), 
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_text("Number of covariates", face = "bold"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))+ 
  # scale_x_reverse()+
  scale_color_brewer(palette = "Set2")





# N fixed ratio vs fdp
ggplot(df_nfixed, aes(x = ratio)) +
  geom_line(aes(y = fdp, color = n), size = 1.5) + 
  geom_hline(yintercept = 0.1, linetype = 2) + 
  labs(y = 'FDP', x = "Ratio p/n") +
  coord_cartesian(ylim = c(0, 0.2)) +
  theme_minimal() + 
  guides(color = guide_legend(title = "Number of observations (Fixed)"))+
  theme(legend.position = c(0.22,0.8), 
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_text("Number of covariates", face = "bold"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))+ 
  scale_color_brewer(palette = "Set2")

# P fixed, ratio vs power
ggplot(df_pfixed, aes(x = ratio)) +
  geom_line(aes(y = power, color = np), size = 1.5) + 
  labs(y = 'Power', x = "Ratio p/n") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() + 
  guides(color = guide_legend(title = "Number of covariates (Fixed)"))+
  theme(legend.position = c(0.21,0.15), 
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_text("Number of covariates", face = "bold"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))+ 
  scale_color_brewer(palette = "Set2")


# n fixed ratio vs power
ggplot(df_nfixed, aes(x = ratio)) +
  geom_line(aes(y = power, color = n), size = 1.5) + 
  labs(y = 'Power', x = "Ratio p/n") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() + 
  guides(color = guide_legend(title = "Number of observations (Fixed)"))+
  theme(legend.position = c(0.76,0.15), 
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_text("Number of covariates", face = "bold"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))+ 
  scale_color_brewer(palette = "Set2")





# P fixed, ratio vs nsel
ggplot(df_pfixed, aes(x = ratio)) +
  geom_line(aes(y = nsel, color = np), size = 1.5) + 
  geom_line(aes(y = k, color = np), linetype = 2) +  
  # geom_hline(yintercept = 0.1, linetype = 2) + 
  labs(y = 'Number of Selected Variables', x = "Ratio p/n") +
  # coord_cartesian(ylim = c(0, 0.2)) +
  theme_minimal() + 
  guides(color = guide_legend(title = "Number of covariates (Fixed)"))+
  theme(legend.position = c(0.2,0.45), 
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_text("Number of covariates", face = "bold"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))+ 
  scale_color_brewer(palette = "Set2")


# N fixed ratio vs nsel
ggplot(df_nfixed, aes(x = ratio)) +
  geom_line(aes(y = nsel, color = n), size = 1.5) + 
  geom_line(aes(y = k, color = n), linetype = 2) + 
  # geom_hline(yintercept = 0.1, linetype = 2) + 
  labs(y = 'Number of Selected Variables', x = "Ratio p/n") +
  coord_cartesian(ylim = c(0, 100)) +
  theme_minimal() + 
  guides(color = guide_legend(title = "Number of observations (Fixed)"))+
  theme(legend.position = c(0.22,0.8), 
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_text("Number of covariates", face = "bold"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))+ 
  scale_color_brewer(palette = "Set2")


##############################







ggplot(df_pfixed, aes(x = ratio)) + 
  geom_line(aes(y = fdp, color = np)) + 
  coord_cartesian(ylim = c(0, 0.12)) 

ggplot(df_nfixed, aes(x = ratio)) + 
  geom_line(aes(y = fdp, color = n)) +
  coord_cartesian(ylim = c(0, 0.12)) 




ggplot(df_nfixed, aes(x = ratio)) +
  geom_line(aes(y = fdp, color = n), size = 1.5) + 
  geom_hline(yintercept = 0.1, linetype = 2) + 
  labs(y = 'FDP', x = expression(paste(rho))) +
  coord_cartesian(ylim = c(0, 0.2)) +
  theme_minimal() + 
  theme(legend.position = c(0.15,0.8), 
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))+ 
  scale_color_brewer(palette = "Set2")


###################################


ggplot(df_pfixed, aes(x = ratio)) + 
  geom_line(aes(y = power, color = np)) + 
  coord_cartesian(ylim = c(0, 1)) 

ggplot(df_pfixed, aes(x = ratio)) + 
  geom_line(aes(y = nsel, color = np)) 


#################3


ggplot(df_nfixed, aes(x = ratio)) + 
  geom_line(aes(y = power, color = n)) +
  coord_cartesian(ylim = c(0, 1)) 
  
ggplot(df_nfixed, aes(x = ratio)) + 
  geom_line(aes(y = nsel, color = n)) 

  # ggplot(sim_1_gaus_plotting, aes(x = corr)) +
  # geom_line(aes(y = nsel, color = Method), size = 1.5) + 
  # labs(y = 'Number of Selected Variables', x = expression(paste(rho))) +
  # theme_minimal() + 
  # theme(legend.position = c(0.15, 0.5), 
  #       legend.background = element_blank(),
  #       legend.box.background = element_rect(colour = "black"),
  #       legend.title = element_text(face = "bold", family = 'Helvetica'),
  #       axis.text.x = element_text(size = 18, family = "Helvetica"),
  #       axis.text.y = element_text(size = 18, family = "Helvetica"),
  #       axis.title.x = element_text(size = 18, family = "Helvetica"),
  #       axis.title.y = element_text(size = 18, family = "Helvetica"))+ 
  # scale_color_brewer(palette = "Set2")

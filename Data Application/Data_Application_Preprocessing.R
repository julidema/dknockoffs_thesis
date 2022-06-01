# ------------------------------------------------
# Real Data Application of Knockoffs
# ------------------------------------------------

library(tidyr)
library(stats)
library(knockoff)
library(dplyr)

# The dataset contains 
# Observations: 61 normal tissue samples and 529 breast cancer samples (Total: 590)
# Covariates: 17,814 genes
# Source: https://data.mendeley.com/datasets/v3cc2p38hb/1

# Read Data
normal <- read.csv("BC-TCGA-Normal.csv", row.names = 1)
tumor <- read.csv("BC-TCGA-Tumor.csv", row.names = 1)

# Check for missing values
sum(is.na(normal))
sum(is.na(tumor))

# No missing values but upon inspection, datasets contain some null values
nnull <- which(normal == 'null', arr.ind = TRUE)
tnull <- which(tumor == 'null', arr.ind = TRUE)

sum(normal == 'null')
sum(tumor == 'null')

nrm <- unique(rownames(nnull))
trm <- unique(rownames(tnull))
uniqueall <- unique(c(nrm, trm))
length(uniqueall) # 534 genes to remove

# Remove these genes containing null values from the datasets
normal_rm <- normal[!(rownames(normal)%in% uniqueall),]
tumor_rm <- tumor[!(rownames(tumor)%in% uniqueall),]


# Create Class category (Normal vs Tumor) (as a row)
normal_rm[dim(normal_rm)[1] + 1, ] <- 0
tumor_rm[dim(tumor_rm)[1] + 1,] <- 1
rownames(normal_rm)[length(rownames(normal_rm))] <- 'Class' 
rownames(tumor_rm)[length(rownames(tumor_rm))] <- 'Class' 


# Transposing the data so that the columns are the variables
# Save names
colnames <- rownames(tumor_rm) # Genes
rownames_t <- colnames(tumor_rm)
rownames_n <- colnames(normal_rm)
rownames <- c(rownames_n, rownames_t) # Subjects

rownames(normal_rm) <- NULL
colnames(normal_rm) <- NULL
rownames(tumor_rm) <- NULL
colnames(tumor_rm) <- NULL


# Combine datasets and transpose 
normal_t <- as.data.frame(t(normal_rm))
tumor_t <- as.data.frame(t(tumor_rm))
all_data_t <- rbind(normal_t, tumor_t)
final_data <- as.data.frame(sapply(all_data_t, as.numeric)) 
rownames(final_data) <- rownames
colnames(final_data) <- colnames


# Standardize data (without the class variable)
std_data <- scale(final_data[,-17281])
data <- cbind(std_data, Class = final_data[,17281])


# Write data to a csv
write.csv(data, "DataApp.csv")

# ------------------------------------------------
# THEN USE THE PYTHON SCRIPT TO PERFORM PCA
# ------------------------------------------------

# Read PCA-reduced data
pcadat <- read.csv("PCA_DataApp.csv", row.names = 1, header = T)

# Read training class
train_class <-read.csv("trainclass.csv", header = TRUE, row.names = 1)

X <- as.matrix(pcadat)
resp <- train_class$Class

# Perform knockoff selection on the PCA reduced dataset
set.seed(1)
X_k <- create.second_order(X, method = 'asdp')
W <- stat.glmnet_coefdiff(X, X_k, resp, family = 'binomial')
t <- knockoff.threshold(W, fdr = 0.1, offset = 0)
selected <- which(W >= t)
selected


# ------------------------------------------------
# THEN USE PYTHON SCRIPT TO PERFORM CLASSIFICATION
# ------------------------------------------------

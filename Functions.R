# ----------------------------------------
# FUNCTIONS USED FOR THESIS SIMULATIONS
# ----------------------------------------
library(knockoff)
library(glmnet)
library(glmnetUtils)

# -------------------------------------------------------------------------------
# Function that mean for every (nonoverlapping) n elements
# -------------------------------------------------------------------------------

comp_mean_n <- function(vec, n){
  test <- c()
  vec[is.na(vec)] <- 0
  
  for(i in 1:(length(vec)/n)){
    test[i] <- mean(vec[(n*(i-1) + 1):(i*n)])  
  }
  return(test)
}

# -------------------------------------------------------------------------------
# Function that computes the number of selected variables
# -------------------------------------------------------------------------------

# ind - column containing index of selected variable (1st column in sel files)
comp_nr_selected <- function(df, ind){
  helper <- diff(df[,ind]) 
  selected <- ifelse(helper[which(helper != 1)] < 0, abs(helper[which(helper != 1)]) + 1, 0)
  last_selection <- ifelse(df[nrow(df),ind+1] == "NONE", 0, df[nrow(df),ind])
  nselected <- c(selected, last_selection)
  return(nselected)
}

# -----------------------------------------------------------------------------------------------
# Function that generates a binary, poisson or gaussian response and applies the knockoff method
# -----------------------------------------------------------------------------------------------

apply_knockoffs <- function(s, n, p, k, a, rho, niter, fdp,
                            family = c("binomial", "poisson", "gaussian"), file1, file2, file3, file4){
  
  calc_fdp <- function(selected) sum(beta[selected] == 0) / max(1, length(selected))

  family = match.arg(family)
  # Data Generation
  if(family == "binomial"){
    set.seed(s)
    mu <- rep(0,p)
    calc_sigma_toep <- function(r, p) {toeplitz(r^(0:(p-1)))}
    sigma_toep <- calc_sigma_toep(rho, p)
    
    calc_X <- function(sigma, p, n) {matrix(rnorm(p*n), n) %*% chol(sigma)}
    X <- calc_X(sigma_toep, p, n)
    
    nonzero <- sample(p, k) 
    beta <- sample(c(-1,1),p, replace = T) * a * (1:p %in% nonzero) / sqrt(n)
    invlogit <- function(x) exp(x) / (1+exp(x))
    y.sample <- function(x) rbinom(n, prob=invlogit(matrix(x, nrow = n) %*% beta), size=1) 
    resp <- factor(y.sample(X), levels = c(0,1), labels = c("A", "B"))
    
    write.table(nonzero, file4, append = TRUE, col.names = FALSE, row.names = TRUE)
    a1 <- data.frame(rho = rho, seed = s, n = n, p = p, k = k, a = a)
    write.table(a1, file1, append = TRUE, col.names = TRUE, row.names = FALSE)
  }
  
  if(family == 'poisson'){
    set.seed(s)
    mu <- rep(0,p)
    calc_sigma_toep <- function(rho, p) {toeplitz(rho^(0:(p-1)))}
    sigma_toep <- calc_sigma_toep(rho, p)
    
    calc_X <- function(sigma, p, n) {matrix(rnorm(p*n), n) %*% chol(sigma)}
    X <- calc_X(sigma_toep, p, n)
    
    nonzero <- sample(p, k) 
    beta <- sample(c(-1,1),p, replace = T) * a * (1:p %in% nonzero) / sqrt(n)
    invlog <- function(x) exp(x) 
    y.sample <- function(x) rpois(n, lambda = invlog(matrix(x, nrow = n) %*% beta)) 
    resp <- y.sample(X)
    
    write.table(nonzero, file4, append = TRUE, col.names = FALSE, row.names = FALSE)
    a1 <- data.frame(rho = rho, seed = s, n = n, p = p, k = k, a = a)
    write.table(a1, file1, append = TRUE, col.names = FALSE)
  }
  
  if(family == 'gaussian'){
    
    set.seed(s)
    mu <- rep(0,p)
    calc_sigma_toep <- function(rho, p) {toeplitz(rho^(0:(p-1)))}
    sigma_toep <- calc_sigma_toep(rho, p)
    
    calc_X <- function(sigma, p, n) {matrix(rnorm(p*n), n) %*% chol(sigma)}
    X <- calc_X(sigma_toep, p, n)
    
    nonzero <- sample(p, k) 
    beta <- sample(c(-1,1),p, replace = T) * a * (1:p %in% nonzero) / sqrt(n)
    y.sample = function(X) X %*% beta + rnorm(n)
    resp <- y.sample(X)
    
    write.table(nonzero, file4, append = TRUE, col.names = FALSE, row.names = FALSE)
    a1 <- data.frame(rho = rho, seed = s, n = n, p = p, k = k, a = a)
    write.table(a1, file1, append = TRUE, col.names = FALSE)
  }
  # Knockoff Procedure
  res <- sapply(1:niter, function(it){
    
    diag_s <- create.solve_sdp(sigma_toep)
    X_k <- create.gaussian(X, mu, sigma_toep, diag_s = diag_s)
    W <- stat.glmnet_coefdiff(X, X_k, resp, family = family)
    t <- knockoff.threshold(W, fdr = fdp, offset = 1)
    selected <- which(W >= t)
    now <- Sys.time()
    a <- data.frame(calc_fdp(selected), sum(beta[selected] != 0)/k, now, fix.empty.names = FALSE)
    b <- data.frame(selected)
    if(nrow(b) == 0){write.table("NONE", file3, append = TRUE, col.names = FALSE)}
    write.table(a, file2, append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(b, file3, append = TRUE, col.names = FALSE, row.names = TRUE)
  })
}

# ---------------------------------------------------------------------------------------------
# Function that generates a binary, Poisson or Gaussian response and applies LASSO with CV
# ---------------------------------------------------------------------------------------------

apply_lasso <- function(s, n, p, k, a, rho, niter,
                        family = c("binomial", "poisson", "gaussian"), file1, file2, file3, file4){
  
  fdp <- function(selected) sum(beta[selected] == 0) / max(1, length(selected))
  
  family = match.arg(family)
  
  # Data Generation  
  if(family == "binomial"){
    set.seed(s)
    mu <- rep(0,p)
    calc_sigma_toep <- function(r, p) {toeplitz(r^(0:(p-1)))}
    sigma_toep <- calc_sigma_toep(rho, p)
    
    calc_X <- function(sigma, p, n) {matrix(rnorm(p*n), n) %*% chol(sigma)}
    X <- calc_X(sigma_toep, p, n)
    
    nonzero <- sample(p, k) # variables that will have nonzero coefficients
    beta <- sample(c(-1,1),p, replace = T) * a * (1:p %in% nonzero) / sqrt(n)
    invlogit <- function(x) exp(x) / (1+exp(x))
    y.sample <- function(x) rbinom(n, prob=invlogit(matrix(x, nrow = n) %*% beta), size=1) # logit link function
    resp <- factor(y.sample(X), levels = c(0,1), labels = c("A", "B"))
    
    write.table(nonzero, file4, append = TRUE, col.names = FALSE, row.names = TRUE)
    a1 <- data.frame(rho = rho, seed = s, n = n, p = p, k = k, a = a)
    write.table(a1, file1, append = TRUE, col.names = TRUE, row.names = FALSE)
  }
  
  if(family == 'poisson'){
    set.seed(s)
    mu <- rep(0,p)
    calc_sigma_toep <- function(rho, p) {toeplitz(rho^(0:(p-1)))}
    sigma_toep <- calc_sigma_toep(rho, p)
    
    calc_X <- function(sigma, p, n) {matrix(rnorm(p*n), n) %*% chol(sigma)}
    X <- calc_X(sigma_toep, p, n)
    
    nonzero <- sample(p, k) # variables that will have nonzero coefficients
    beta <- sample(c(-1,1),p, replace = T) * a * (1:p %in% nonzero) / sqrt(n)
    invlog <- function(x) exp(x) 
    y.sample <- function(x) rpois(n, lambda = invlog(matrix(x, nrow = n) %*% beta)) 
    resp <- y.sample(X)
    
    write.table(nonzero, file4, append = TRUE, col.names = FALSE, row.names = FALSE)
    a1 <- data.frame(rho = rho, seed = s, n = n, p = p, k = k, a = a)
    write.table(a1, file1, append = TRUE, col.names = FALSE)
  }
  
  if(family == 'gaussian'){
    
    set.seed(s)
    mu <- rep(0,p)
    calc_sigma_toep <- function(rho, p) {toeplitz(rho^(0:(p-1)))}
    sigma_toep <- calc_sigma_toep(rho, p)
    
    calc_X <- function(sigma, p, n) {matrix(rnorm(p*n), n) %*% chol(sigma)}
    X <- calc_X(sigma_toep, p, n)
    
    nonzero <- sample(p, k) # variables that will have nonzero coefficients
    beta <- sample(c(-1,1),p, replace = T) * a * (1:p %in% nonzero) / sqrt(n)
    y.sample = function(X) X %*% beta + rnorm(n)
    resp <- y.sample(X)
    
    fdp <- function(selected) sum(beta[selected] == 0) / max(1, length(selected))
    
    write.table(nonzero, file4, append = TRUE, col.names = FALSE, row.names = FALSE)
    a1 <- data.frame(rho = rho, seed = s, n = n, p = p, k = k, a = a)
    write.table(a1, file1, append = TRUE, col.names = FALSE)
  }
  
  
  # Apply Lasso
  res <- sapply(1:niter, function(it){
    
    lasso_cv <- cv.glmnet(X, resp, family = family, parallel = T) # 10-fold cross validation
    lambda_best <- lasso_cv$lambda.min
    lasso_best <- glmnet(X, resp, family = family, lambda = lambda_best)
    result <- coef(lasso_best) 
    selected <- which(result != 0)[-1] - 1 # Exclude intercept and shift by 1
    now <- Sys.time()
    a <- data.frame(fdp(selected), sum(beta[selected] != 0)/k, now, fix.empty.names = FALSE)
    b <- data.frame(selected)
    if(nrow(b) == 0){write.table("NONE", file3, append = TRUE, col.names = FALSE)}
    write.table(a, file2, append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(b, file3, append = TRUE, col.names = FALSE, row.names = TRUE)
    
  })
}

# ---------------------------------------------------------------------------------------------
# Function that generates a binary, Poisson or Gaussian response and applies Elastic Net with CV
# ---------------------------------------------------------------------------------------------

apply_elastic_net <- function(s, n, p, k, a, rho, niter,
                              family = c("binomial", "poisson", "gaussian"), file1, file2, file3, file4){
  
  fdp <- function(selected) sum(beta[selected] == 0) / max(1, length(selected))
  
  family = match.arg(family)
  
  # Data Generation  
  if(family == "binomial"){
    set.seed(s)
    mu <- rep(0,p)
    calc_sigma_toep <- function(r, p) {toeplitz(r^(0:(p-1)))}
    sigma_toep <- calc_sigma_toep(rho, p)
    
    calc_X <- function(sigma, p, n) {matrix(rnorm(p*n), n) %*% chol(sigma)}
    X <- calc_X(sigma_toep, p, n)
    
    nonzero <- sample(p, k) # variables that will have nonzero coefficients
    beta <- sample(c(-1,1),p, replace = T) * a * (1:p %in% nonzero) / sqrt(n)
    invlogit <- function(x) exp(x) / (1+exp(x))
    y.sample <- function(x) rbinom(n, prob=invlogit(matrix(x, nrow = n) %*% beta), size=1) # logit link function
    resp <- factor(y.sample(X), levels = c(0,1), labels = c("A", "B"))
    
    write.table(nonzero, file4, append = TRUE, col.names = FALSE, row.names = TRUE)
    a1 <- data.frame(rho = rho, seed = s, n = n, p = p, k = k, a = a)
    write.table(a1, file1, append = TRUE, col.names = TRUE, row.names = FALSE)
  }
  
  if(family == 'poisson'){
    set.seed(s)
    mu <- rep(0,p)
    calc_sigma_toep <- function(rho, p) {toeplitz(rho^(0:(p-1)))}
    sigma_toep <- calc_sigma_toep(rho, p)
    
    calc_X <- function(sigma, p, n) {matrix(rnorm(p*n), n) %*% chol(sigma)}
    X <- calc_X(sigma_toep, p, n)
    
    nonzero <- sample(p, k) # variables that will have nonzero coefficients
    beta <- sample(c(-1,1),p, replace = T) * a * (1:p %in% nonzero) / sqrt(n)
    invlog <- function(x) exp(x) 
    y.sample <- function(x) rpois(n, lambda = invlog(matrix(x, nrow = n) %*% beta)) 
    resp <- y.sample(X)
    
    write.table(nonzero, file4, append = TRUE, col.names = FALSE, row.names = FALSE)
    a1 <- data.frame(rho = rho, seed = s, n = n, p = p, k = k, a = a)
    write.table(a1, file1, append = TRUE, col.names = FALSE)
  }
  
  if(family == 'gaussian'){
    
    set.seed(s)
    mu <- rep(0,p)
    calc_sigma_toep <- function(rho, p) {toeplitz(rho^(0:(p-1)))}
    sigma_toep <- calc_sigma_toep(rho, p)
    
    calc_X <- function(sigma, p, n) {matrix(rnorm(p*n), n) %*% chol(sigma)}
    X <- calc_X(sigma_toep, p, n)
    
    nonzero <- sample(p, k) # variables that will have nonzero coefficients
    beta <- sample(c(-1,1),p, replace = T) * a * (1:p %in% nonzero) / sqrt(n)
    y.sample = function(X) X %*% beta + rnorm(n)
    resp <- y.sample(X)
    
    fdp <- function(selected) sum(beta[selected] == 0) / max(1, length(selected))
    
    write.table(nonzero, file4, append = TRUE, col.names = FALSE, row.names = FALSE)
    a1 <- data.frame(rho = rho, seed = s, n = n, p = p, k = k, a = a)
    write.table(a1, file1, append = TRUE, col.names = FALSE)
  }
  
  
  # Apply Elastic Net with alpha = 0.5 (Fixed) and cross validation for lambda
  res <- sapply(1:niter, function(it){
    
    elastic_cv <- cv.glmnet(X, resp, family = family, parallel = T, alpha = 0.5) # 10-fold cross validation for lambda
    lambda_best <- elastic_cv$lambda.min
    elastic_best <- glmnet(X, resp, family = family, lambda = lambda_best, alpha = 0.5)
    result <- coef(elastic_best) 
    selected <- which(result != 0)[-1] - 1 # Exclude intercept and shift by 1 to account for the intercept
    now <- Sys.time()
    a <- data.frame(fdp(selected), sum(beta[selected] != 0)/k, now, fix.empty.names = FALSE)
    b <- data.frame(selected)
    if(nrow(b) == 0){write.table("NONE", file3, append = TRUE, col.names = FALSE)}
    write.table(a, file2, append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(b, file3, append = TRUE, col.names = FALSE, row.names = TRUE)
    
  })
}



# ---------------------------------------------------------------------------------------------
# Function that generates a binary, Poisson or Gaussian response and applies adaptive Lasso
# ---------------------------------------------------------------------------------------------

apply_lasso_adaptive <- function(s, n, p, k, a, rho, niter,
                                 family = c("binomial", "poisson", "gaussian"), file1, file2, file3, file4){
  
  fdp <- function(selected) sum(beta[selected] == 0) / max(1, length(selected))
  family = match.arg(family)
  
  # Data Generation  
  if(family == "binomial"){
    set.seed(s)
    mu <- rep(0,p)
    calc_sigma_toep <- function(r, p) {toeplitz(r^(0:(p-1)))}
    sigma_toep <- calc_sigma_toep(rho, p)
    
    calc_X <- function(sigma, p, n) {matrix(rnorm(p*n), n) %*% chol(sigma)}
    X <- calc_X(sigma_toep, p, n)
    
    nonzero <- sample(p, k) # variables that will have nonzero coefficients
    beta <- sample(c(-1,1),p, replace = T) * a * (1:p %in% nonzero) / sqrt(n)
    invlogit <- function(x) exp(x) / (1+exp(x))
    y.sample <- function(x) rbinom(n, prob=invlogit(matrix(x, nrow = n) %*% beta), size=1) # logit link function
    resp <- factor(y.sample(X), levels = c(0,1), labels = c("A", "B"))
    
    data_for_weights <- data.frame(resp, X)
    glm_out <- glm(resp ~., data_for_weights, family = 'binomial')
    weights <- 1/abs(coef(glm_out))[-1]
    
    write.table(nonzero, file4, append = TRUE, col.names = FALSE, row.names = TRUE)
    a1 <- data.frame(rho = rho, seed = s, n = n, p = p, k = k, a = a)
    write.table(a1, file1, append = TRUE, col.names = TRUE, row.names = FALSE)
  }
  
  if(family == 'poisson'){
    set.seed(s)
    mu <- rep(0,p)
    calc_sigma_toep <- function(rho, p) {toeplitz(rho^(0:(p-1)))}
    sigma_toep <- calc_sigma_toep(rho, p)
    
    calc_X <- function(sigma, p, n) {matrix(rnorm(p*n), n) %*% chol(sigma)}
    X <- calc_X(sigma_toep, p, n)
    
    nonzero <- sample(p, k) # variables that will have nonzero coefficients
    beta <- sample(c(-1,1),p, replace = T) * a * (1:p %in% nonzero) / sqrt(n)
    invlog <- function(x) exp(x) 
    y.sample <- function(x) rpois(n, lambda = invlog(matrix(x, nrow = n) %*% beta)) 
    resp <- y.sample(X)
    
    data_for_weights <- data.frame(resp, X)
    glm_out <- glm(resp ~., data_for_weights, family = 'poisson')
    weights <- 1/abs(coef(glm_out))[-1]
    
    write.table(nonzero, file4, append = TRUE, col.names = FALSE, row.names = FALSE)
    a1 <- data.frame(rho = rho, seed = s, n = n, p = p, k = k, a = a)
    write.table(a1, file1, append = TRUE, col.names = FALSE)
  }
  
  if(family == 'gaussian'){
    
    set.seed(s)
    mu <- rep(0,p)
    calc_sigma_toep <- function(rho, p) {toeplitz(rho^(0:(p-1)))}
    sigma_toep <- calc_sigma_toep(rho, p)
    
    calc_X <- function(sigma, p, n) {matrix(rnorm(p*n), n) %*% chol(sigma)}
    X <- calc_X(sigma_toep, p, n)
    
    nonzero <- sample(p, k) # variables that will have nonzero coefficients
    beta <- sample(c(-1,1),p, replace = T) * a * (1:p %in% nonzero) / sqrt(n)
    y.sample = function(X) X %*% beta + rnorm(n)
    resp <- y.sample(X)
    
    # Weights for adaptive lasso
    data_for_weights <- data.frame(resp, X)
    lm_out <- lm(resp ~., data_for_weights)
    weights <- 1/abs(coef(lm_out))[-1]
    
    write.table(nonzero, file4, append = TRUE, col.names = FALSE, row.names = FALSE)
    a1 <- data.frame(rho = rho, seed = s, n = n, p = p, k = k, a = a)
    write.table(a1, file1, append = TRUE, col.names = FALSE)
  }
  
  
  res <- sapply(1:niter, function(it){
    
    lasso_cv <- cv.glmnet(X, resp, family = family, parallel = T, penalty.factor = weights) # 10-fold cross validation
    lambda_best <- lasso_cv$lambda.min
    lasso_best <- glmnet(X, resp, family = family, lambda = lambda_best)
    result <- coef(lasso_best) 
    selected <- which(result != 0)[-1] - 1 # Exclude intercept and shift by 1
    now <- Sys.time()
    a <- data.frame(fdp(selected), sum(beta[selected] != 0)/k, now, fix.empty.names = FALSE)
    b <- data.frame(c("START",selected))
    if(nrow(b) == 0){write.table("NONE", file3, append = TRUE, col.names = FALSE)}
    write.table(a, file2, append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(b, file3, append = TRUE, col.names = FALSE, row.names = TRUE)
    
  })
}



# ----------------------------------------------------------------------------------------------
# Function that generates a Binary, Poisson or Gaussian response and applies simultaneously 
# both Knockoffs and Knockoffs+ at various levels of q
# ----------------------------------------------------------------------------------------------
apply_knockoffs_orp <- function(s, n, p, k, a, rho, niter, 
                                family = c("binomial", "poisson", "gaussian"), file1, file2, file3, file4,
                                file5, file6, file7, file8, file9, file10, file11, file12, file13, file14,
                                file15, file16, file17, file18, file19, file20, file21, file22){
  
  calc_fdp <- function(selected) sum(beta[selected] == 0) / max(1, length(selected))
  
  family = match.arg(family)
  if(family == "binomial"){
    set.seed(s)
    mu <- rep(0,p)
    calc_sigma_toep <- function(r, p) {toeplitz(r^(0:(p-1)))}
    sigma_toep <- calc_sigma_toep(rho, p)
    
    calc_X <- function(sigma, p, n) {matrix(rnorm(p*n), n) %*% chol(sigma)}
    X <- calc_X(sigma_toep, p, n)
    
    nonzero <- sample(p, k) # variables that will have nonzero coefficients
    beta <- sample(c(-1,1),p, replace = T) * a * (1:p %in% nonzero) / sqrt(n)
    invlogit <- function(x) exp(x) / (1+exp(x))
    y.sample <- function(x) rbinom(n, prob=invlogit(matrix(x, nrow = n) %*% beta), size=1) # logit link function
    resp <- factor(y.sample(X), levels = c(0,1), labels = c("A", "B"))
    
    write.table(nonzero, file22, append = TRUE, col.names = FALSE, row.names = TRUE)
    a1 <- data.frame(rho = rho, seed = s, n = n, p = p, k = k, a = a)
    write.table(a1, file1, append = TRUE, col.names = TRUE, row.names = FALSE)
  }
  
  if(family == 'poisson'){
    set.seed(s)
    mu <- rep(0,p)
    calc_sigma_toep <- function(rho, p) {toeplitz(rho^(0:(p-1)))}
    sigma_toep <- calc_sigma_toep(rho, p)
    
    calc_X <- function(sigma, p, n) {matrix(rnorm(p*n), n) %*% chol(sigma)}
    X <- calc_X(sigma_toep, p, n)
    
    nonzero <- sample(p, k) # variables that will have nonzero coefficients
    beta <- sample(c(-1,1),p, replace = T) * a * (1:p %in% nonzero) / sqrt(n)
    invlog <- function(x) exp(x) 
    y.sample <- function(x) rpois(n, lambda = invlog(matrix(x, nrow = n) %*% beta)) 
    resp <- y.sample(X)
    
    write.table(nonzero, file22, append = TRUE, col.names = FALSE, row.names = FALSE)
    a1 <- data.frame(rho = rho, seed = s, n = n, p = p, k = k, a = a)
    write.table(a1, file1, append = TRUE, col.names = FALSE)
  }
  
  if(family == 'gaussian'){
    
    set.seed(s)
    mu <- rep(0,p)
    calc_sigma_toep <- function(rho, p) {toeplitz(rho^(0:(p-1)))}
    sigma_toep <- calc_sigma_toep(rho, p)
    
    calc_X <- function(sigma, p, n) {matrix(rnorm(p*n), n) %*% chol(sigma)}
    X <- calc_X(sigma_toep, p, n)
    
    nonzero <- sample(p, k) # variables that will have nonzero coefficients
    beta <- sample(c(-1,1),p, replace = T) * a * (1:p %in% nonzero) / sqrt(n)
    y.sample = function(X) X %*% beta + rnorm(n)
    resp <- y.sample(X)
    
    write.table(nonzero, file22, append = TRUE, col.names = FALSE, row.names = FALSE)
    a1 <- data.frame(rho = rho, seed = s, n = n, p = p, k = k, a = a)
    write.table(a1, file1, append = TRUE, col.names = FALSE)
  }
  
  res <- sapply(1:niter, function(it){
    
    diag_s <- create.solve_sdp(sigma_toep)
    X_k <- create.gaussian(X, mu, sigma_toep, diag_s = diag_s)
    W <- stat.glmnet_coefdiff(X, X_k, resp, family = family)
    
    # Knockoffs+ thresholds
    t_plus1 <- knockoff.threshold(W, fdr = 0.05, offset = 1)
    t_plus2 <- knockoff.threshold(W, fdr = 0.1, offset = 1)
    t_plus3 <- knockoff.threshold(W, fdr = 0.2, offset = 1)
    t_plus4 <- knockoff.threshold(W, fdr = 0.3, offset = 1)
    t_plus5 <- knockoff.threshold(W, fdr = 0.5, offset = 1)
    # Knockoffs thresholds
    t1 <- knockoff.threshold(W,fdr = 0.05, offset = 0)
    t2 <- knockoff.threshold(W,fdr = 0.1, offset = 0)
    t3 <- knockoff.threshold(W,fdr = 0.2, offset = 0)
    t4 <- knockoff.threshold(W,fdr = 0.3, offset = 0)
    t5 <- knockoff.threshold(W,fdr = 0.5, offset = 0)
    
    # Knockoffs+ selected
    selected_plus1 <- which(W >= t_plus1)
    selected_plus2 <- which(W >= t_plus2)
    selected_plus3 <- which(W >= t_plus3)
    selected_plus4 <- which(W >= t_plus4)
    selected_plus5 <- which(W >= t_plus5)
    # Knockoffs selected
    selected1 <- which(W >= t1)
    selected2 <- which(W >= t2)
    selected3 <- which(W >= t3)
    selected4 <- which(W >= t4)
    selected5 <- which(W >= t5)
    
  
    now <- Sys.time()
  
    # Write results to files
    a_plus1 <- data.frame(calc_fdp(selected_plus1), sum(beta[selected_plus1] != 0)/k, now, fix.empty.names = FALSE)
    b_plus1 <- data.frame(selected_plus1)
    a_plus2 <- data.frame(calc_fdp(selected_plus2), sum(beta[selected_plus2] != 0)/k, now, fix.empty.names = FALSE)
    b_plus2 <- data.frame(selected_plus2)
    a_plus3 <- data.frame(calc_fdp(selected_plus3), sum(beta[selected_plus3] != 0)/k, now, fix.empty.names = FALSE)
    b_plus3 <- data.frame(selected_plus3)
    a_plus4 <- data.frame(calc_fdp(selected_plus4), sum(beta[selected_plus4] != 0)/k, now, fix.empty.names = FALSE)
    b_plus4 <- data.frame(selected_plus4)
    a_plus5 <- data.frame(calc_fdp(selected_plus5), sum(beta[selected_plus5] != 0)/k, now, fix.empty.names = FALSE)
    b_plus5 <- data.frame(selected_plus5)
    
    a1 <- data.frame(calc_fdp(selected1), sum(beta[selected1] != 0)/k, now, fix.empty.names = FALSE)
    b1 <- data.frame(selected1)
    a2 <- data.frame(calc_fdp(selected2), sum(beta[selected2] != 0)/k, now, fix.empty.names = FALSE)
    b2 <- data.frame(selected2)
    a3 <- data.frame(calc_fdp(selected3), sum(beta[selected3] != 0)/k, now, fix.empty.names = FALSE)
    b3 <- data.frame(selected3)
    a4 <- data.frame(calc_fdp(selected4), sum(beta[selected4] != 0)/k, now, fix.empty.names = FALSE)
    b4 <- data.frame(selected4)
    a5 <- data.frame(calc_fdp(selected5), sum(beta[selected5] != 0)/k, now, fix.empty.names = FALSE)
    b5 <- data.frame(selected5)
    
    
    if(nrow(b1) == 0){write.table("NONE", file3, append = TRUE, col.names = FALSE)}
    write.table(a1, file2, append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(b1, file3, append = TRUE, col.names = FALSE, row.names = TRUE)
    
    if(nrow(b2) == 0){write.table("NONE", file5, append = TRUE, col.names = FALSE)}
    write.table(a2, file4, append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(b2, file5, append = TRUE, col.names = FALSE, row.names = TRUE)
    
    if(nrow(b3) == 0){write.table("NONE", file7, append = TRUE, col.names = FALSE)}
    write.table(a3, file6, append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(b3, file7, append = TRUE, col.names = FALSE, row.names = TRUE)
    
    if(nrow(b4) == 0){write.table("NONE", file9, append = TRUE, col.names = FALSE)}
    write.table(a4, file8, append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(b4, file9, append = TRUE, col.names = FALSE, row.names = TRUE)
    
    if(nrow(b5) == 0){write.table("NONE", file11, append = TRUE, col.names = FALSE)}
    write.table(a5, file10, append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(b5, file11, append = TRUE, col.names = FALSE, row.names = TRUE)
    
    
    if(nrow(b_plus1) == 0){write.table("NONE", file13, append = TRUE, col.names = FALSE)}
    write.table(a_plus1, file12, append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(b_plus1, file13, append = TRUE, col.names = FALSE, row.names = TRUE)
    
    if(nrow(b_plus2) == 0){write.table("NONE", file15, append = TRUE, col.names = FALSE)}
    write.table(a_plus2, file14, append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(b_plus2, file15, append = TRUE, col.names = FALSE, row.names = TRUE)
    
    if(nrow(b_plus3) == 0){write.table("NONE", file17, append = TRUE, col.names = FALSE)}
    write.table(a_plus3, file16, append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(b_plus3, file17, append = TRUE, col.names = FALSE, row.names = TRUE)
    
    if(nrow(b_plus4) == 0){write.table("NONE", file19, append = TRUE, col.names = FALSE)}
    write.table(a_plus4, file18, append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(b_plus4, file19, append = TRUE, col.names = FALSE, row.names = TRUE)
    
    if(nrow(b_plus5) == 0){write.table("NONE", file21, append = TRUE, col.names = FALSE)}
    write.table(a_plus5, file20, append = TRUE, col.names = FALSE, row.names = FALSE)
    write.table(b_plus5, file21, append = TRUE, col.names = FALSE, row.names = TRUE)
    
  })
}

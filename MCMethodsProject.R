
#########################
# GENERAL CONFIGURATION #
#########################

# LIBRARIES
library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)
library(scales)
library(patchwork)

# Common parameters
nsim <- 10000           # No. of simulations for an empirical characterization
n <- 1:500              # Vector of sample sizes
gamma <- 3              # Indicator for Pr(X > gamma)

##########################
# RAW MONTE CARLO N(0,1) #
##########################

# Simulation of target distribution N(0,1)
# Estimation of E[X] and Pr(X > gamma)
# Characterization: bias, variance y MSE

# Inicialization and result structure
results <- data.frame(n = integer(), 
                      estimator = character(), 
                      value = numeric(),
                      bias = numeric(), 
                      variance = numeric(), 
                      mse = numeric())
results_samples <- list(Mean = list(), Probs = list())

# Simmulation and characterization for every samplesize
pb <- txtProgressBar(min = 0, max = length(n), style = 3)
set.seed(1234)
for (j in 1:length(n)) {
  
  mean_estimates <- numeric(nsim)
  prob_estimates <- numeric(nsim)
  
  for (i in 1:nsim) {
    samples <- rnorm(n = n[j], mean = 0, sd = 1)  # N(0, 1) samples
    # Mean estimation
    mean_estimates[i] <- mean(samples)
    # Estimation of P(X > gamma)
    prob_estimates[i] <- mean(samples > gamma)
  }
  
  # Theoretic Values
  true_mean <- 0
  true_prob <- pnorm(gamma, lower.tail = FALSE)
  
  # Characterizing estimator
  # Mean
  bias_mean <- mean(mean_estimates) - true_mean
  var_mean <- var(mean_estimates)
  mse_mean <- bias_mean^2 + var_mean
  # P(X > gamma)
  bias_prob <- mean(prob_estimates) - true_prob
  var_prob <- var(prob_estimates)
  mse_prob <- bias_prob^2 + var_prob
  
  # Results
  results[j, ] <- data.frame(n = n[j], estimator = "Mean", 
                             value = mean(mean_estimates), bias = bias_mean, 
                             variance = var_mean, mse = mse_mean)
  results[j + length(n), ] <- data.frame(n = n[j], estimator = "Prob",
                                         value = mean(prob_estimates),
                                         bias = bias_prob, variance = var_prob,
                                         mse = mse_prob)
  results_samples$Mean[[j]] <- mean_estimates
  results_samples$Probs[[j]] <- prob_estimates
  setTxtProgressBar(pb, j)
}
close(pb)

############################
# PLOTS OF RAW MONTE CARLO #
############################

# Histogram of estimators for N = 20
# Bias, variance, and MSE vs. sample size
# Estimation trajectory over n

# p20_mean, p20_prob: distributions of estimates
# p1, p2: convergence of estimates over n
# pbias_*, pvar_*, pmse_*: performance metrics

p20_mean <- ggplot(as_tibble(results_samples$Mean[[20]]), aes(x = value)) +
  geom_histogram(alpha = 0.7, fill = "gray95", color = "black", size = 0.1) +
  geom_vline(xintercept = 0, color = "black", size = 0.1) +
  geom_vline(xintercept = mean(results_samples$Mean[[20]]), color = "black", 
             size = 0.1, linetype = "dashed") +
  labs(x = "Estimate", y = "Frequency", title = expression("E[" * X * "]")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("label", x = 0.485, y = 990, 
           label = paste("Bias:", sprintf("%.3f", results[results$n == 20 & results$estimator == "Mean", "bias"]),
                         "\nVariance:", sprintf("%.3f", results[results$n == 20 & results$estimator == "Mean", "variance"]),
                         "\nMSE:", sprintf("%.3f", results[results$n == 20 & results$estimator == "Mean", "mse"])),
           size = 3, hjust = 0)

p20_prob <- ggplot(as_tibble(results_samples$Probs[[20]]), aes(x = value)) +
  geom_histogram(alpha = 0.7, fill = "gray95", color = "black", size = 0.1) +
  geom_vline(xintercept = pnorm(gamma, lower.tail = FALSE), color = "black", 
             size = 0.1) +
  geom_vline(xintercept = mean(results_samples$Probs[[20]]), color = "black", 
             size = 0.1, linetype = "dashed") +
  labs(x = "Estimate", y = "Frequency", title = expression(Pr(X > gamma))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("label", x = 0.115, y = 8900,
           label = paste("Bias:", sprintf("%.4f", results[results$n == 20 & results$estimator == "Prob", "bias"]),
                         "\nVariance:", sprintf("%.4f", results[results$n == 20 & results$estimator == "Prob", "variance"]),
                         "\nMSE:", sprintf("%.4f", results[results$n == 20 & results$estimator == "Prob", "mse"])),
           size = 3, hjust = 0)

p20_mean + p20_prob + plot_layout(ncol = 2)

# Sample size vs Estimator
p1 <- ggplot(as_tibble(results[results$estimator == "Mean", c("n", "value")]), 
             aes(x = n, y = value)) +
  geom_line(size = 0.1) +
  geom_hline(yintercept = 0, color = "red") +
  labs(x = "Sample size", y = "Estimator") +
  ggtitle("E(X)") +
  theme_bw()

p2 <- ggplot(as_tibble(results[results$estimator == "Prob", c("n", "value")]), 
             aes(x = n, y = value)) +
  geom_line(size = 0.1) +
  geom_hline(yintercept = pnorm(gamma, lower.tail = FALSE), col = "red") +
  labs(x = "Sample size", y = "Estimator") +
  ggtitle("P(X > γ)") +
  theme_bw()

p1 + p2 + plot_layout(nrow = 1)

# Bias, variance and MSE
pbias_mean <- ggplot(results[results$estimator == "Mean", ], 
                     aes(x = n, y = bias)) +
  geom_line(size = 0.1) +
  scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  labs(x = "Sample size", y = "Bias", title = expression("E[" * X * "]")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))

pvar_mean <- ggplot(results[results$estimator == "Mean", ], 
                    aes(x = n, y = variance)) +
  geom_line(size = 0.1) +
  scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  labs(x = "Sample size", y = "Variance") +
  theme_bw()

pmse_mean <- ggplot(results[results$estimator == "Mean", ], 
                    aes(x = n, y = mse)) +
  geom_line(size = 0.1) +
  scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  labs(x = "Sample size", y = "MSE") +
  theme_bw()

pbias_prob <- ggplot(results[results$estimator == "Prob", ], 
                     aes(x = n, y = bias)) +
  geom_line(size = 0.1) +
  scale_y_continuous(labels = label_number(accuracy = 0.0001)) +
  labs(x = "Sample size", y = "Bias", title = expression(Pr(X>gamma))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))

pvar_prob <- ggplot(results[results$estimator == "Prob", ], 
                    aes(x = n, y = variance)) +
  geom_line(size = 0.1) +
  scale_y_continuous(labels = label_number(accuracy = 0.0001)) +
  labs(x = "Sample size", y = "Variance") +
  theme_bw()

pmse_prob <- ggplot(results[results$estimator == "Prob", ], 
                    aes(x = n, y = mse)) +
  geom_line(size = 0.1) +
  scale_y_continuous(labels = label_number(accuracy = 0.0001)) +
  labs(x = "Sample size", y = "MSE") +
  theme_bw()

pbias_mean + pbias_prob + pvar_mean + pvar_prob + pmse_mean + pmse_prob +
  plot_layout(nrow = 3, ncol = 2, guides = "collect") &
  theme(legend.position = "right",         
        legend.justification = c(1, 1),       
        legend.box.just = "top")

################################
# IMPORTANCE SAMPLING OPTION 2 #
################################

# Target: 0.5 N(-3, 1) + 0.5 N(5, 4)
# Good proposal: N(0, 25)
# Bad proposal:  N(-3, 25)

# For each proposal:
#   - sample x ~ q
#   - compute p(x), q(x), and weights
#   - compute IS estimate of E[X]
#   - visualize weights and density overlap

#################
# GOOD PROPOSAL #
#################

# sample x ~ q
x <- rnorm(n = 1000, mean = 0, sd = sqrt(25))

# compute p(x) and q(x)
p_x <- 0.5 * dnorm(x = x, mean = -3, sd = 1) +    # Target
  0.5 * dnorm(x = x, mean = 5, sd = 4)
q_x <- dnorm(x = x, mean = 0, sd = sqrt(25))      # Proposal

# Weights
w <- p_x / q_x

# Normalized weights
w_normalized <- w / sum(w)

# IS Estimator
is <- sum(x * w_normalized)

# Plot
gp <- ggplot(data.frame(x = x, w_scaled = w_normalized * max(q_x) / 
                    quantile(w_normalized, 0.735)), 
       aes(x = x, y = w_scaled)) +
  geom_segment(aes(xend = x, y = 0, yend = w_scaled), color = "lightgray", 
               alpha = 0.7) +
  geom_point(x = x, y = 0, size = 0.5) +
  stat_function(fun = function(x) {
    0.5 * dnorm(x, mean = -3, sd = 1) +
      0.5 * dnorm(x, mean = 5, sd = 4)
  }, color = "steelblue") +
  stat_function(fun = dnorm, args = list(mean = 0, sd = sqrt(25)),
                color = "darkorange", linetype = "dashed") +
  xlim(-20, 20) +
  labs(y = "Density", title = "N(0, 25)") +
  theme_bw()

################
# BAD PROPOSAL #
################

# sample x ~ q
x <- rnorm(n = 1000, mean = -3, sd = sqrt(25))

# compute p(x) and q(x)
p_x <- 0.5 * dnorm(x = x, mean = -3, sd = 1) +    # Target
  0.5 * dnorm(x = x, mean = 5, sd = 4)
q_x <- dnorm(x = x, mean = -3, sd = sqrt(25))     # Proposal

# Weights
w <- p_x / q_x

# Normalized weights
w_normalized <- w / sum(w)

# IS Estimator
is <- sum(x * w_normalized)

# Plot
bp <- ggplot(data.frame(x = x, w_scaled = w_normalized * max(q_x) / 
                    quantile(w_normalized, 0.65)), 
       aes(x = x, y = w_scaled)) +
  geom_segment(aes(xend = x, y = 0, yend = w_scaled), color = "lightgray", 
               alpha = 0.7) +
  geom_point(x = x, y = 0, size = 0.5) +
  stat_function(fun = function(x) {
    0.5 * dnorm(x, mean = -3, sd = 1) +
      0.5 * dnorm(x, mean = 5, sd = 4)
  }, color = "steelblue") +
  stat_function(fun = dnorm, args = list(mean = -3, sd = sqrt(25)),
                color = "darkorange", linetype = "dashed") +
  xlim(-22, 20) +
  labs(y = "Density", title = "N(-3, 25)") +
  theme_bw()

gp + bp

#####################################################################
#  IMPORTANCE SAMPLING WITH BIMODAL TARGET – PERFORMANCE EVALUATION #
#####################################################################

# Target: 0.5 N(-3,1) + 0.5 N(5,4)
# Proposal: N(0, 25)

# Estimators are:
#   - E[X], E[X^2], Pr(X > gamma)
#   - UIS and SNIS
# Evaluate over different sample sizes
# Empirical bias, variance, MSE

# Common parameters
nsim <- 10000           # No. of simulations for an empirical characterization
n <- 1:500              # Vector of sample sizes
gamma <- 2              # Indicator for Pr(X > gamma)

results <- data.frame(n = integer(), 
                      estimator = character(), 
                      type = character(),
                      value = numeric(),
                      bias = numeric(), 
                      variance = numeric(), 
                      mse = numeric())
results_samples <- list(Mean = list(), SecMom = list(), Probs = list())

# Simmulation and characterization for every samplesize
pb <- txtProgressBar(min = 0, max = length(n), style = 3)
set.seed(1234)
for (j in 1:length(n)) {
  
  mean_estimates <- matrix(NA, nrow = nsim, ncol = 2)
  secmom_estimates <- matrix(NA, nrow = nsim, ncol = 2)
  prob_estimates <- matrix(NA, nrow = nsim, ncol = 2)
  
  for (i in 1:nsim) {
    
    # Samples of q(x)
    x <- rnorm(n = n[j], mean = 0, sd = sqrt(25))
    
    # Densities
    p_x <- 0.5 * dnorm(x = x, mean = -3, sd = 1) +    # Target
      0.5 * dnorm(x = x, mean = 5, sd = 4)
    q_x <- dnorm(x = x, mean = 0, sd = sqrt(25))      # Proposal
    
    # Weights
    w <- p_x / q_x
    
    # Normalized weights
    w_normalized <- w / sum(w)
    
    # Mean Estimator
    mean_estimates[i, 1:2] <- c(mean(x * w), sum(x * w_normalized))
    
    # Second moment Estimator
    secmom_estimates[i, 1:2] <- c(mean(x^2 * w), sum(x^2 * w_normalized))
    
    # Estimator of P(X > gamma)
    prob_estimates[i, 1:2] <- c(mean(as.numeric(x > gamma) * w),
                                sum(as.numeric(x > gamma) * w_normalized))
  }
  
  # Theoretical values
  true_mean <- 1
  true_secmom <- 25.5
  true_prob <- pnorm(gamma, lower.tail = FALSE)
  
  # Estimator characterization
  # Mean
  bias_mean <- apply(mean_estimates, 2, mean) - true_mean
  var_mean <- apply(mean_estimates, 2, var)
  mse_mean <- bias_mean^2 + var_mean
  
  # Second moment
  bias_secmom <- apply(secmom_estimates, 2, mean) - true_secmom
  var_secmom <- apply(secmom_estimates, 2, var)
  mse_secmom <- bias_secmom^2 + var_secmom
  
  # P(X > gamma)
  bias_prob <- apply(prob_estimates, 2, mean) - true_prob
  var_prob <- apply(prob_estimates, 2, var)
  mse_prob <- bias_prob^2 + var_prob
  
  # Results
  results[(6*(j-1) + 1):(6*(j-1) + 2), ] <- data.frame(n = n[j], 
                                                       estimator = "Mean", 
                                                       type = c("UIS", "SNIS"),
                                                       value = apply(mean_estimates, 2, mean),
                                                       bias = bias_mean, 
                                                       variance = var_mean, 
                                                       mse = mse_mean)
  results[(6*(j-1) + 3):(6*(j-1) + 4), ] <- data.frame(n = n[j], 
                                                       estimator = "SecMom",
                                                       type = c("UIS", "SNIS"),
                                                       value = apply(secmom_estimates, 2, mean),
                                                       bias = bias_secmom, 
                                                       variance = var_secmom,
                                                       mse = mse_secmom)
  results[(6*(j-1) + 5):(6*(j-1) + 6), ] <- data.frame(n = n[j], 
                                                       estimator = "Prob",
                                                       type = c("UIS", "SNIS"),
                                                       value = apply(prob_estimates, 2, mean),
                                                       bias = bias_prob, 
                                                       variance = var_prob,
                                                       mse = mse_prob)
  results_samples$Mean[[j]] <- mean_estimates
  results_samples$SecMom[[j]] <- secmom_estimates
  results_samples$Probs[[j]] <- prob_estimates
  setTxtProgressBar(pb, j)
}
close(pb)

results$estimator <- factor(results$estimator,
                            levels = c("Mean", "SecMom", "Prob"))
results$type <- factor(results$type,
                       levels = c("UIS", "SNIS"))

###########################################
# PLOTS OF IS MONTE CARLO BINOMIAL TARGET #
###########################################

# Histogram of estimators for N = 20
# Bias, variance, and MSE vs. sample size
# Estimation trajectory over n

# p20_mean, p20_prob: distributions of estimates
# p1, p2: convergence of estimates over n
# pbias_*, pvar_*, pmse_*: performance metrics

pbias_mean <- ggplot(results[results$estimator == "Mean", ], 
                     aes(x = n, y = bias, linetype = type)) +
  geom_line(size = 0.1) +
  scale_color_brewer(
    palette = "Set1",
    labels = c(Mean = expression(E * "[" * X * "]"),
               SecMom = expression(E * "[" * X^2 * "]"),
               Prob = expression(Pr(X>gamma)))) +
  labs(x = "Sample size", y = "Bias", color = "Parameter", 
       linetype = "Estimator", title = expression(E * "[" * X * "]")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))

pbias_secmom <- ggplot(results[results$estimator == "SecMom", ], 
                       aes(x = n, y = bias, linetype = type)) +
  geom_line(size = 0.1) +
  scale_color_brewer(
    palette = "Set1",
    labels = c(Mean = expression(E * "[" * X * "]"),
               SecMom = expression(E * "[" * X^2 * "]"),
               Prob = expression(Pr(X>gamma)))) +
  labs(x = "Sample size", y = "Bias", color = "Parameter", 
       linetype = "Estimator", title = expression(E * "[" * X^2 * "]")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))

pbias_prob <- ggplot(results[results$estimator == "Prob", ], 
                     aes(x = n, y = bias, linetype = type)) +
  geom_line(size = 0.1) +
  scale_color_brewer(
    palette = "Set1",
    labels = c(Mean = expression(E * "[" * X * "]"),
               SecMom = expression(E * "[" * X^2 * "]"),
               Prob = expression(Pr(X>gamma)))) +
  labs(x = "Sample size", y = "Bias", color = "Parameter", 
       linetype = "Estimator", title = expression(Pr(X>gamma))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))

pvar_mean <- ggplot(results[results$estimator == "Mean", ], 
                    aes(x = n, y = variance, linetype = type)) +
  geom_line(size = 0.1) +
  scale_color_brewer(
    palette = "Set1",
    labels = c(Mean = expression(E * "[" * X * "]"),
               SecMom = expression(E * "[" * X^2 * "]"),
               Prob = expression(Pr(X>gamma)))) +
  labs(x = "Sample size", y = "Variance", color = "Parameter", 
       linetype = "Estimator") +
  theme_bw() +
  theme(legend.position = "none")

pvar_secmom <- ggplot(results[results$estimator == "SecMom", ], 
                      aes(x = n, y = variance, linetype = type)) +
  geom_line(size = 0.1) +
  scale_color_brewer(
    palette = "Set1",
    labels = c(Mean = expression(E * "[" * X * "]"),
               SecMom = expression(E * "[" * X^2 * "]"),
               Prob = expression(Pr(X>gamma)))) +
  labs(x = "Sample size", y = "Variance", color = "Parameter", 
       linetype = "Estimator") +
  theme_bw() +
  theme(legend.position = "none")

pvar_prob <- ggplot(results[results$estimator == "Prob", ], 
                    aes(x = n, y = variance, linetype = type)) +
  geom_line(size = 0.1) +
  scale_color_brewer(
    palette = "Set1",
    labels = c(Mean = expression(E * "[" * X * "]"),
               SecMom = expression(E * "[" * X^2 * "]"),
               Prob = expression(Pr(X>gamma)))) +
  labs(x = "Sample size", y = "Variance", color = "Parameter", 
       linetype = "Estimator") +
  theme_bw() +
  theme(legend.position = "none")

pmse_mean <- ggplot(results[results$estimator == "Mean", ], 
                    aes(x = n, y = mse, linetype = type)) +
  geom_line(size = 0.1) +
  scale_color_brewer(
    palette = "Set1",
    labels = c(Mean = expression(E * "[" * X * "]"),
               SecMom = expression(E * "[" * X^2 * "]"),
               Prob = expression(Pr(X>gamma)))) +
  labs(x = "Sample size", y = "MSE", color = "Parameter", 
       linetype = "Estimator") +
  theme_bw() +
  theme(legend.position = "none")

pmse_secmom <- ggplot(results[results$estimator == "SecMom", ], 
                      aes(x = n, y = mse, linetype = type)) +
  geom_line(size = 0.1) +
  scale_color_brewer(
    palette = "Set1",
    labels = c(Mean = expression(E * "[" * X * "]"),
               SecMom = expression(E * "[" * X^2 * "]"),
               Prob = expression(Pr(X>gamma)))) +
  labs(x = "Sample size", y = "MSE", color = "Parameter", 
       linetype = "Estimator") +
  theme_bw() +
  theme(legend.position = "none")

pmse_prob <- ggplot(results[results$estimator == "Prob", ], 
                    aes(x = n, y = mse, linetype = type)) +
  geom_line(size = 0.1) +
  scale_color_brewer(
    palette = "Set1",
    labels = c(Mean = expression(E * "[" * X * "]"),
               SecMom = expression(E * "[" * X^2 * "]"),
               Prob = expression(Pr(X>gamma)))) +
  labs(x = "Sample size", y = "MSE", color = "Parameter", 
       linetype = "Estimator") +
  theme_bw() +
  theme(legend.position = "none")

pbias_mean + pbias_secmom + pbias_prob +
  pvar_mean + pvar_secmom + pvar_prob +
  pmse_mean + pmse_secmom + pmse_prob +
  plot_layout(nrow = 3, ncol = 3, guides = "collect") &
  theme(legend.position = "right",         
        legend.justification = c(1, 1),       
        legend.box.just = "top")

#############################################
# SENSITIVITY TO PROPOSAL MEAN AND VARIANCE #
#############################################

# Parameters
mu_grid <- c(-3, 1, 4)
sigma_grid <- c(2, sqrt(24.5), 8)
n <- 20
nsim <- 10000
gamma <- 2

# Ground truth for bimodal target: mixture of N(-3,1) and N(5,4)
true_mean <- 1
true_secmom <- 25.5
true_prob <- pnorm(gamma, lower.tail = FALSE)

# Store all results
results_all <- list()

# Loop over proposal combinations
for (mu_q in mu_grid) {
  for (sigma_q in sigma_grid) {
    
    results <- data.frame()
    pb <- txtProgressBar(min = 0, max = length(n), style = 3)
    set.seed(1234)
    
    for (j in seq_along(n)) {
      mean_estimates <- matrix(NA, nrow = nsim, ncol = 2)
      secmom_estimates <- matrix(NA, nrow = nsim, ncol = 2)
      prob_estimates <- matrix(NA, nrow = nsim, ncol = 2)
      
      for (i in 1:nsim) {
        x <- rnorm(n = n[j], mean = mu_q, sd = sigma_q)
        
        p_x <- 0.5 * dnorm(x, mean = -3, sd = 1) + 0.5 * dnorm(x, mean = 5, 
                                                               sd = 4)
        q_x <- dnorm(x, mean = mu_q, sd = sigma_q)
        w <- p_x / q_x
        w_normalized <- w / sum(w)
        
        mean_estimates[i, ] <- c(mean(x * w), sum(x * w_normalized))
        secmom_estimates[i, ] <- c(mean(x^2 * w), sum(x^2 * w_normalized))
        prob_estimates[i, ] <- c(mean((x > gamma) * w), 
                                 sum((x > gamma) * w_normalized))
      }
      
      # Compute statistics
      bias_mean <- colMeans(mean_estimates) - true_mean
      var_mean <- apply(mean_estimates, 2, var)
      mse_mean <- bias_mean^2 + var_mean
      
      bias_secmom <- colMeans(secmom_estimates) - true_secmom
      var_secmom <- apply(secmom_estimates, 2, var)
      mse_secmom <- bias_secmom^2 + var_secmom
      
      bias_prob <- colMeans(prob_estimates) - true_prob
      var_prob <- apply(prob_estimates, 2, var)
      mse_prob <- bias_prob^2 + var_prob
      
      # Combine results
      results <- rbind(
        results,
        data.frame(n = n[j], estimator = "Mean", type = c("UIS", "SNIS"),
                   bias = bias_mean, variance = var_mean, mse = mse_mean),
        data.frame(n = n[j], estimator = "SecMom", type = c("UIS", "SNIS"),
                   bias = bias_secmom, variance = var_secmom, mse = mse_secmom),
        data.frame(n = n[j], estimator = "Prob", type = c("UIS", "SNIS"),
                   bias = bias_prob, variance = var_prob, mse = mse_prob)
      )
      setTxtProgressBar(pb, j)
    }
    close(pb)
    
    results$proposal <- paste0("mu=", mu_q, "_sigma=", sigma_q)
    results$estimator <- factor(results$estimator, 
                                levels = c("Mean", "SecMom", "Prob"))
    results$type <- factor(results$type, levels = c("UIS", "SNIS"))
    
    results_all[[paste0("mu=", mu_q, "_sigma=", sigma_q)]] <- results
  }
}

# Combine all proposal results
results_combined <- bind_rows(results_all)

library(dplyr)
library(flextable)
library(officer)

# Set filter parameters
target_type <- "SNIS"
target_n <- n

# Prepare summary table
summary_all <- results_combined |>
  filter(n == target_n, type == target_type) |>
  separate(proposal, into = c("mu", "sigma"), sep = "_") |>
  mutate(
    mu = as.numeric(gsub("mu=", "", mu)),
    sigma = as.numeric(gsub("sigma=", "", sigma))
  ) |>
  rename(
    `Proposal Mean` = mu,
    `Proposal SD` = sigma,
    Estimator = estimator,
    Bias = bias,
    Variance = variance,
    MSE = mse
  ) |>
  arrange(Estimator, `Proposal Mean`, `Proposal SD`) |>
  select(Estimator, `Proposal Mean`, `Proposal SD`, Bias, Variance, MSE)

highlighted_table <- summary_all |>
  group_by(Estimator) |>
  mutate(
    highlight_bias = abs(Bias) == min(abs(Bias), na.rm = TRUE),
    highlight_var  = abs(Variance) == min(abs(Variance), na.rm = TRUE),
    highlight_mse  = abs(MSE) == min(abs(MSE), na.rm = TRUE)
  ) |>
  ungroup()

summary_all <- summary_all |>
  mutate(
    Bias = round(Bias, 3),
    Variance = round(Variance, 3),
    MSE = round(MSE, 3),
    `Proposal SD` = ceiling(`Proposal SD`)  # round up to nearest whole number
  )

summary_all <- summary_all |>
  mutate(
    Estimator = recode(Estimator,
                       "Mean" = "E(X)",
                       "SecMom" = "E(X^2)",
                       "Prob" = "P(X > γ)"
    )
  )


# Define borders
outer_border <- fp_border(color = "black", width = 2)
inner_border <- fp_border(color = "gray70", width = 1)

# Reuse highlight_table from earlier if it includes `min_bias`, etc.
# Or just use summary_all if you dropped that logic

ft <- flextable(summary_all) |>
  # Math-style estimator labels
  compose(j = "Estimator", value = as_paragraph(as_i(Estimator))) |>
  
  # Basic style
  set_caption(paste0("Performance Metrics for all estimators (", target_type, 
                     ", n = ", target_n, ")")) |>
  autofit() |>
  theme_booktabs() |>
  align(align = "center", part = "all") |>
  colformat_num(j = c("Bias", "Variance", "MSE"), digits = 4) |>
  
  # Font styling
  bold(part = "header") |>
  font(part = "all", fontname = "Inconsolata") |>
  fontsize(size = 10, part = "body") |>
  
  # Grouped rows
  merge_v(j = c("Estimator", "Proposal Mean")) |>
  
  # Borders and grid styling
  border_remove() |>
  border_outer(border = outer_border, part = "all") |>
  border_inner_v(border = inner_border, part = "body") |>
  border_inner_h(border = inner_border, part = "body")

# Define custom red border
red_border <- fp_border(color = "red", width = 1.5)

# Apply borders manually to the cell block at row 5
ft |>
  hline(i = 4, j = 4:6, border = red_border, part = "body") |>   # top edge
  hline(i = 5, j = 4:6, border = red_border, part = "body") |>   # bottom edge
  vline(i = 5, j = 3, border = red_border, part = "body") |>     # left edge
  vline(i = 5, j = 6, border = red_border, part = "body")        # right edge


###########################
# IS IN HIGHER DIMENSIONS #
###########################

# Parameters
d_vec <- seq(1:5)
n <- 20                  # number of samples
nsim <- 10000            # number of simulations per dimension
mu_q <- 1                # mean of proposal
sigma_q <- sqrt(24.5)    # std dev of proposal
gamma <- 2               # threshold for tail probability

# Data frame to store MSE by dimension
mse_by_dim <- data.frame()

# Store a sample of weights for visualization at highest dimension
weights_snapshot <- NULL

pb <- txtProgressBar(min = 0, max = length(d_vec), style = 3)

set.seed(1234)
for (j in seq_along(d_vec)) {
  d <- d_vec[j]
  
  estimates <- matrix(NA, nrow = nsim, ncol = 3)
  colnames(estimates) <- c("E(X)", "E(X^2)", "P(X > γ)")
  
  for (i in 1:nsim) {
    # Sample x dimension by dimension (same as vectorized sampling here due to independence)
    x <- matrix(NA, nrow = n, ncol = d)
    for (k in 1:d) {
      x[, k] <- rnorm(n, mean = mu_q, sd = sigma_q)
    }
    
    # Log densities of target (mixture of Gaussians)
    log_p_x <- rowSums(log(
      0.5 * dnorm(x, mean = -3, sd = 1) + 0.5 * dnorm(x, mean = 5, sd = 4)
    ))
    
    # Log densities of proposal
    log_q_x <- rowSums(dnorm(x, mean = mu_q, sd = sigma_q, log = TRUE))
    
    # Importance weights (unnormalized then normalized)
    log_w <- log_p_x - log_q_x
    w <- exp(log_w)
    w_norm <- w / sum(w)
    
    # Store normalized weights at highest dimension for a snapshot
    if (j == length(d_vec) && i == 1) {
      weights_snapshot <- sort(w_norm, decreasing = TRUE)
    }
    
    # Estimators over the first coordinate
    x1 <- x[, 1]
    estimates[i, 1] <- sum(x1 * w_norm)
    estimates[i, 2] <- sum(x1^2 * w_norm)
    estimates[i, 3] <- sum((x1 > gamma) * w_norm)
  }
  
  # True values for comparison
  true_mean <- 1
  true_secmom <- 25.5
  true_prob <- integrate(
    function(x) (0.5 * dnorm(x, -3, 1) + 0.5 * dnorm(x, 5, 4)) * (x > gamma),
    -Inf, Inf
  )$value
  
  true_vals <- c(true_mean, true_secmom, true_prob)
  mse_vals <- colMeans((estimates - matrix(true_vals, nrow = nsim, ncol = 3, byrow = TRUE))^2)
  
  mse_by_dim <- rbind(
    mse_by_dim,
    data.frame(
      dimension = d,
      parameter = c("E(X)", "E(X^2)", "P(X > γ)"),
      mse = mse_vals
    )
  )
  
  setTxtProgressBar(pb, j)
}

close(pb)

#################################
# Plotting the MSE vs dimension #
#################################

library(ggplot2)
library(dplyr)
library(patchwork)

base_theme <- theme_bw(base_family = "Inconsolata") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )

# E(X)
p_mean <- ggplot(filter(mse_by_dim, parameter == "E(X)"),
                 aes(x = dimension, y = mse)) +
  geom_line(color = "black", size = 1) +
  scale_x_continuous(breaks = d_vec) +
  labs(title = expression(MSE~"for "~E(X)),
       x = "Dimension", y = "MSE") +
  base_theme

# E(X²)
p_secmom <- ggplot(filter(mse_by_dim, parameter == "E(X^2)"),
                   aes(x = dimension, y = mse)) +
  geom_line(color = "black", size = 1) +
  scale_x_continuous(breaks = d_vec) +
  labs(title = expression(MSE~"for "~E(X^2)),
       x = "Dimension", y = "MSE") +
  base_theme

# P(X > γ)
p_prob <- ggplot(filter(mse_by_dim, parameter == "P(X > γ)"),
                 aes(x = dimension, y = mse)) +
  geom_line(color = "black", size = 1) +
  scale_x_continuous(breaks = d_vec) +
  labs(title = expression(MSE~"for "~P(X > gamma)),
       x = "Dimension", y = "MSE") +
  base_theme

# Combined MSE plot
p_mean / p_secmom / p_prob

####################################################
# Plot sorted normalized weights at d = max(d_vec) #
####################################################

df_weights <- data.frame(index = seq_along(weights_snapshot),
                         weight = weights_snapshot)

ggplot(df_weights, aes(x = index, y = weight)) +
  geom_line(color = "black", size = 1) +
  scale_x_continuous(breaks = 1:20) +
  labs(title = paste("Sorted normalized IS weights at d =", max(d_vec)),
       x = "Sample index (ordered by weight)", y = "Normalized weight") +
  base_theme


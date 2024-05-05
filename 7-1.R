# Load necessary library
library(MASS) # for multivariate normal distribution

# Set parameters for the simulation
set.seed(123)
n <- 100 # number of observations
p <- 100  # number of parameters (dimensions)
true_theta <- rep(5, p) # true parameter vector
cov_matrix <- diag(p) # covariance matrix for the distribution

# Simulate data
X <- mvrnorm(n = n, mu = true_theta, Sigma = cov_matrix)

# MLE Estimator
mle_estimates <- colMeans(X)

# James-Stein Estimator
js_shrinkage <- function(estimates, n, p) {
  mean_estimates <- mean(estimates)
  S_squared <- sum((estimates - mean_estimates)^2) / n
  shrinkage_factor <- 1 - ((p - 2) / S_squared)
  shrinkage_factor <- max(min(shrinkage_factor, 1), 0)  # Ensuring the shrinkage factor is between 0 and 1
  shrunk_estimates <- shrinkage_factor * estimates + (1 - shrinkage_factor) * mean_estimates
  return(shrunk_estimates)
}

js_estimates <- js_shrinkage(mle_estimates, n, p)

# Calculate Mean Squared Error
mse_mle <- sum((mle_estimates - true_theta)^2) / p
mse_js <- sum((js_estimates - true_theta)^2) / p

# Print the results
print(paste("MSE for MLE: ", mse_mle))
print(paste("MSE for James-Stein: ", mse_js))

# Compare estimators
if (mse_js < mse_mle) {
  print("James-Stein estimator has lower MSE.")
} else {
  print("MLE has lower MSE.")
}


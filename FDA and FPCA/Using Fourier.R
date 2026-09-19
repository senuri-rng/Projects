# Functional Data Analysis and FPCA
# Project 1: Fourier Basis and Smoothing

# Load prepared data
source("R/01_data_preparation.R")

# Number of Fourier basis functions
M <- 3

# Construct Fourier basis manually
fourier_basis_mat <- function(t, M, period) {
  
  basis <- matrix(1, nrow = length(t), ncol = 1)
  
  for (m in 1:((M - 1) / 2)) {
    basis <- cbind(
      basis,
      sin(2 * pi * m * t / period),
      cos(2 * pi * m * t / period)
    )
  }
  
  basis
}

# Create Fourier basis
B <- fourier_basis_mat(t, M, period = P)

# Second-order difference penalty
D <- diff(diag(P), differences = 2)
R <- t(D %*% B) %*% (D %*% B)

# Smoothing parameters
lambda_grid <- 10^seq(-5, 5, length.out = 20)
gcv_values <- numeric(length(lambda_grid))

# Select smoothing parameter using GCV
for (i in seq_along(lambda_grid)) {
  
  lambda <- lambda_grid[i]
  
  C_temp <- solve(
    t(B) %*% B + lambda * R,
    t(B) %*% as.matrix(Y)
  )
  
  Y_hat_temp <- B %*% C_temp
  
  hat_matrix <- B %*%
    solve(t(B) %*% B + lambda * R, t(B))
  
  df <- sum(diag(hat_matrix))
  
  gcv_values[i] <-
    mean(colSums((Y - Y_hat_temp)^2)) /
    (1 - df / P)^2
}

# Optimal smoothing parameter
best_lambda_fourier <- lambda_grid[which.min(gcv_values)]

# Final penalized least-squares fit
C_hat_fourier <- solve(
  t(B) %*% B + best_lambda_fourier * R,
  t(B) %*% as.matrix(Y)
)

Y_hat_fourier <- B %*% C_hat_fourier

# Display selected smoothing parameter
best_lambda_fourier

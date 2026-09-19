# Fourier Basis and Penalized Smoothing

# Number of Fourier basis functions
M = 5

# Time points
t = 1:P

# Fourier basis function
fourier_basis_mat = function(t, M, period){
  basis = matrix(1, nrow = length(t), ncol = 1)
  
  for(m in 1:((M-1)/2)){
    basis = cbind(
      basis,
      sin(2*pi*m*t / period),
      cos(2*pi*m*t / period)
    )
  }
  
  basis
}

# Fourier basis matrix
B = fourier_basis_mat(t, M, period = P)
dim(B)

# Penalized Smoothing

# Second-order difference penalty
D = diff(diag(P), differences = 2)
R = t(D %*% B) %*% (D %*% B)

# Smoothing parameter grid
lambda_grid = 10^seq(-3, 3, length.out = 20)
gcv_values = numeric(length(lambda_grid))

# Generalized Cross-Validation
for(i in seq_along(lambda_grid)){
  
  lambda = lambda_grid[i]
  
  C_hat = solve( t(B) %*% B + lambda * R, t(B) %*% as.matrix(Y))
  
  Y_hat = B %*% C_hat
  
  hat_matrix = B %*% solve(t(B) %*% B + lambda * R, t(B))
  
  df = sum(diag(hat_matrix))
  
  gcv_values[i] = mean(colSums((Y - Y_hat)^2)) / (1 - df/P)^2
}

# Optimal smoothing parameter
best_lambda = lambda_grid[which.min(gcv_values)]

# Final fitted curves
C_hat = solve( t(B) %*% B + best_lambda * R) %*% t(B) %*% as.matrix(Y)

Y_hat = B %*% C_hat

# Functional Principal Component Analysis (FPCA)

# Numerical integration weights
w = rep(0, P)

w[1:(P-1)] <- w[1:(P-1)] + 0.5 * (t[-1] - t[-P])

w[2:P] <- w[2:P] + 0.5 * (t[-1] - t[-P])

# Weight matrix
W = diag(w)

# Inner-product matrix
J = t(B) %*% W %*% B

# Eigen-Decomposition of Inner-Product Matrix

eigJ = eigen(J, symmetric = TRUE)

valsJ = eigJ$values
vecsJ = eigJ$vectors

Jhalf = vecsJ %*% diag(sqrt(valsJ)) %*% t(vecsJ)

# Covariance Matrix

Mmat = (1/K) * Jhalf %*% (C_hat %*% t(C_hat)) %*% Jhalf

# Functional Principal Components

eigMmat = eigen(Mmat, symmetric = TRUE)

rho = eigMmat$values
U = eigMmat$vectors

JhalfInv = vecsJ %*% diag((valsJ)^(-1/2)) %*% t(vecsJ)

Bcoef = JhalfInv %*% U

# Functional principal components
g = B %*% Bcoef

# Proportion of Variance Explained

pve = rho / sum(rho)
cum_pve = cumsum(pve)

# Functional Principal Component Plots

par(mfrow = c(ceiling(ncol(g)/2), 2))

for (i in 1:ncol(g)) {
  
  plot(
    t,
    g[, i],
    type = "l",
    lwd = 2,
    col = "blue",
    main = paste("Functional Principal Component", i),
    xlab = "Time",
    ylab = "Value"
  )
}

# FPC Scores

scores = t(Bcoef) %*% J %*% C_hat

dim(scores)

# Eigenfunctions and FPC Scores

matplot(
  t,
  g,
  type = "l",
  lty = 1,
  main = "Eigenfunctions (Functional PCs)",
  xlab = "t",
  ylab = "g(t)"
)

matplot(
  t(scores[1:3, ]),
  type = "b",
  main = "FPC Scores (First 3 Components)",
  xlab = "Function Index",
  ylab = "Score Value"
)

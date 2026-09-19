# Fourier Basis and Penalized Smoothing
# Without using fda package

Data = read.csv("/Users/senurigunaratne/Desktop/Research Work/R Codes/Project 1/Canadian_climate_history.csv")
Precipitation_Data = Data[,c(3,5,7,9,11,13,15,17,19,21,23,25,27)]

Y = na.omit(Precipitation_Data) 

#P: no. of rows in Y dataset
#K: no. of cols in Y dataset
#M: no of bases

P = nrow(Y)
K = ncol(Y)
t = 1:P
#t = seq(0, 1, length=P)
M = 5

# Fourier basis function
fourier_basis_mat = function(t, M, period){
  basis = matrix(1, nrow = length(t), ncol = 1)  #1st column is 1
  for(m in 1:((M-1)/2)){
    basis = cbind(basis, sin(2*pi*m*t / period), cos(2*pi*m*t / period))
  }
  basis
}

# Fourier basis matrix
B = fourier_basis_mat(t, M, period = P); B  
dim(B)

# Second-order difference penalty
D = diff(diag(P), differences = 2)
R = t(D%*%B) %*% (D%*%B)  

# Smoothing parameter grid
lambda_grid = 10^seq(-3, 3, length.out = 20)
gcv_values = numeric(length(lambda_grid))

# Generalized Cross-Validation
for(i in seq_along(lambda_grid)){
  lambda = lambda_grid[i]
  C_hat = solve(t(B) %*% B + lambda*R, t(B) %*% as.matrix(Y))
  Y_hat = B %*% C_hat
  hat_matrix = B %*% solve(t(B) %*% B + lambda*R, t(B))
  df = sum(diag(hat_matrix))
  gcv_values[i] = mean(colSums((Y - Y_hat)^2)) / (1 - df/P)^2
}

# Optimal smoothing parameter
best_lambda = lambda_grid[which.min(gcv_values)]; best_lambda

# Final fitted curves
C_hat = solve(t(B)%*%B + best_lambda*R)%*% t(B) %*% as.matrix(Y)
Y_hat = B %*% C_hat 

# Functional Principal Component Analysis (FPCA)
# Numerical integration weights
#0.5*(t[-1]-t[-P])

w = rep(0, P)
w[1:(P-1)] <- w[1:(P-1)] + 0.5*(t[-1]-t[-P])
w[2:P] <- w[2:P] + 0.5*(t[-1]-t[-P])

#w = rep(1, P)
#w[1] = 0.5 
#w[P] = 0.5

# Weight matrix
W = diag(w)
dim(W)

# Inner-product matrix
J = t(B)%*%W%*%B

# Eigen-decomposition of J
eigJ = eigen(J, symmetric = TRUE)
valsJ = eigJ$values
vecsJ = eigJ$vectors

# Matrix square root of J
Jhalf = vecsJ%*%diag(sqrt(valsJ))%*%t(vecsJ)

# Covariance matrix
Mmat = (1/K)*Jhalf%*%(C_hat%*%t(C_hat))%*%Jhalf

# Eigen-decomposition of covariance matrix
eigMmat = eigen(Mmat, symmetric = TRUE)
rho = eigMmat$values
U = eigMmat$vectors

# Inverse square root of J
JhalfInv = vecsJ%*%diag((valsJ)^(-1/2))%*%t(vecsJ)

# Functional principal component coefficients
Bcoef = JhalfInv%*%U

#t(Bcoef)%*%J%*%Bcoef #should be equal to I
#t(U)%*%U #should be equal to I
#dim(B%*%Bcoef)

# Functional principal components
g = B%*%Bcoef

#plot(1:nrow(g), g[,3], type="l")
#plot(g)
#matplot(g)

# Proportion of Variance Explained
pve = rho / sum(rho)
cum_pve = cumsum(pve)

# Plot Functional Principal Components
#png("pca_fourier.png", width = 1000, height = 400)

par(mfrow = c(ceiling(ncol(g)/2), 2))  # adjust layout

for (i in 1:ncol(g)) {
  plot(t, g[, i], type = "l", lwd = 2, col = "blue",
       main = paste("Functional Principal Component", i),
       xlab = "Time", ylab = "Value")
}

#dev.off() 

# FPC Scores
scores = t(Bcoef) %*% J %*% C_hat 
dim(scores)

# Plot Eigenfunctions and FPC Scores
matplot(t, g, type = "l", lty = 1, main = "Eigenfunctions (Functional PCs)",
        xlab = "t", ylab = "g(t)")

matplot(t(scores[1:3, ]), type = "b", main = "FPC Scores (first 3 components)",
        xlab = "Function index", ylab = "Score value")

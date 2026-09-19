# Data Preparation
Data = read.csv("/Users/senurigunaratne/Desktop/Research Work/R Codes/Project 1/Canadian_climate_history.csv")
Precipitation_Data = Data[, c(3,5,7,9,11,13,15,17,19,21,23,25,27)]
Y = na.omit(Precipitation_Data)

P = nrow(Y)
K = ncol(Y)
t = seq(0, 1, length.out = P)

M = 10
degree = 3 # cubic spline, order = 4

# B-Spline Basis Using splines Package
library(splines)

B = bs(t, df = M, degree = degree, intercept = TRUE)
dim(B)

# Penalized Smoothing Using splines Package
D = diff(diag(P), differences = 2)
R = t(D%*%B) %*% (D%*%B)  # penalty matrix

lambda_grid = 10^seq(-5, 5, length.out = 20)
gcv_values = numeric(length(lambda_grid))

for(i in seq_along(lambda_grid)){
  lambda = lambda_grid[i]
  C_temp = solve(t(B) %*% B + lambda*R, t(B) %*% as.matrix(Y))
  Y_hat_temp = B %*% C_temp
  hat_matrix = B %*% solve(t(B) %*% B + lambda*R, t(B))
  df = sum(diag(hat_matrix))
  gcv_values[i] = mean(colSums((Y - Y_hat_temp)^2)) / (1 - df/P)^2
}

best_lambda = lambda_grid[which.min(gcv_values)]; best_lambda

C_hat = solve(t(B)%*%B + best_lambda*R)%*% t(B) %*% as.matrix(Y)
Y_hat = B %*% C_hat

# Manual B-Spline Basis Using Cox-de Boor
in_knots = M - degree - 1
internal_knots = seq(0, 1, length.out = in_knots + 2)[2:(in_knots + 1)]
knots = c(rep(0, degree + 1), internal_knots, rep(1, degree + 1))


# Using Cox-de Boor
# t = time points where you want to evaluate the spline
# k = index of the basis function
# d = degree of the spline
# knots = the knot sequence (boundary knots repeated degree+1 times)
b_spline_basis_recursive <- function(t, k, d, knots) {
  if (d == 0) { # constant function
    left <- knots[k]
    right <- knots[k + 1]
    return(ifelse((t >= left & t < right) | (right == max(knots) & t == right), 1, 0))
  } else {
    denom1 <- knots[k + d] - knots[k]
    denom2 <- knots[k + d + 1] - knots[k + 1]
    
    term1 <- if (denom1 == 0) rep(0, length(t)) else ((t - knots[k]) / denom1) * b_spline_basis_recursive(t, k, d - 1, knots)
    term2 <- if (denom2 == 0) rep(0, length(t)) else ((knots[k + d + 1] - t) / denom2) * b_spline_basis_recursive(t, k + 1, d - 1, knots)
    
    return(term1 + term2)
  }
}

n_bases = length(knots) - degree - 1
B_rec = sapply(1:n_bases, function(i) b_spline_basis_recursive(t, i, degree, knots))
dim(B_rec)

# Comparison of B-Spline Basis Functions
# Plotting
#png("bsplines_comparison_2.png", width = 1000, height = 400)
par(mfrow = c(1, 2))

matplot(t, B, type = "l", lty = 1, col = "blue",
        main = "B-spline Basis splines package",
        xlab = "Time", ylab = "Basis Values")

matplot(t, B_rec, type = "l", lty = 1, col = "red",
        main = "B-spline Basis constructed manually",
        xlab = "Time", ylab = "Basis Values")

#dev.off()
               
# Penalized Smoothing Using Manual B-Spline Basis
D = diff(diag(P), differences = 2)
R = t(D%*%B_rec) %*% (D%*%B_rec)  # penalty matrix

lambda_grid = 10^seq(-5, 5, length.out = 20)
gcv_values = numeric(length(lambda_grid))

for(i in seq_along(lambda_grid)){
  lambda = lambda_grid[i]
  C_temp = solve(t(B_rec) %*% B_rec + lambda*R, t(B_rec) %*% as.matrix(Y))
  Y_hat_temp = B_rec %*% C_temp
  hat_matrix = B_rec %*% solve(t(B_rec) %*% B_rec + lambda*R, t(B_rec))
  df = sum(diag(hat_matrix))
  gcv_values[i] = mean(colSums((Y - Y_hat_temp)^2)) / (1 - df/P)^2
}

best_lambda_Spline = lambda_grid[which.min(gcv_values)]; best_lambda_Spline

C_hat_spline = solve(t(B_rec)%*%B_rec + best_lambda_Spline*R)%*% t(B_rec) %*% as.matrix(Y)
Y_hat_spline = B_rec %*% C_hat_spline

# Comparison of Smoothing Results
col_index <- 1

#png("bsplines_comparison.png", width = 1000, height = 400)

par(mfrow = c(1, 2))

# Left plot: manual smoothing result
plot(t, Y_hat[, col_index],
     type = "l", col = "blue", lwd = 2,
     xlab = "Time", ylab = "Precipitation",
     main = "Manual BSpline Smoothing")

# Right plot: FDA package result
plot(t, Y_hat_spline[, col_index],
     type = "l", col = "red", lwd = 2, lty = 2,
     xlab = "Time", ylab = "Precipitation",
     main = "Splines Package Smoothing")

#dev.off()
               
# FPCA Using splines Package
# Data preparation for FPCA

Data = read.csv("/Users/senurigunaratne/Desktop/Research Work/Fall 2025/R code_Projects/Project 1/Canadian_climate_history.csv")
Precipitation_Data = Data[, c(3,5,7,9,11,13,15,17,19,21,23,25,27)]
Y = na.omit(Precipitation_Data)

P = nrow(Y)
K = ncol(Y)
t = seq(0, 1, length.out = P)

M = 10
degree = 3 # cubic spline, order = 4


# B-Spline Basis for FPCA
library(splines)

B = bs(t, df = M, degree = degree, intercept = TRUE)
dim(B)

D = diff(diag(P), differences = 2)
R = t(D%*%B) %*% (D%*%B)  # penalty matrix

lambda_grid = 10^seq(-5, 5, length.out = 20)
gcv_values = numeric(length(lambda_grid))

for(i in seq_along(lambda_grid)){
  lambda = lambda_grid[i]
  C_temp = solve(t(B) %*% B + lambda*R, t(B) %*% as.matrix(Y))
  Y_hat_temp = B %*% C_temp
  hat_matrix = B %*% solve(t(B) %*% B + lambda*R, t(B))
  df = sum(diag(hat_matrix))
  gcv_values[i] = mean(colSums((Y - Y_hat_temp)^2)) / (1 - df/P)^2
}

best_lambda = lambda_grid[which.min(gcv_values)]; best_lambda

C_hat = solve(t(B)%*%B + best_lambda*R)%*% t(B) %*% as.matrix(Y)
Y_hat = B %*% C_hat

# Numerical Integration and Inner Product Matrix
#0.5*(t[-1]-t[-P])

w = rep(0, P)
w[1:(P-1)] <- w[1:(P-1)] + 0.5*(t[-1]-t[-P])
w[2:P] <- w[2:P] + 0.5*(t[-1]-t[-P])

#w = rep(1, P)
#w[1] = 0.5 
#w[P] = 0.5

W = diag(w)
dim(W)

J = t(B)%*%W%*%B

# Eigen-Decomposition of Inner Product Matrix
eigJ = eigen(J, symmetric = TRUE)
valsJ = eigJ$values
vecsJ = eigJ$vectors

Jhalf = vecsJ%*%diag(sqrt(valsJ))%*%t(vecsJ)
               
# Covariance Matrix and FPCA
Mmat = (1/K)*Jhalf%*%(C_hat%*%t(C_hat))%*%Jhalf

eigMmat = eigen(Mmat, symmetric = TRUE)
rho = eigMmat$values
U = eigMmat$vectors

JhalfInv = vecsJ%*%diag((valsJ)^(-1/2))%*%t(vecsJ)

Bcoef = JhalfInv%*%U

#t(Bcoef)%*%J%*%Bcoef #should be equal to I
#t(U)%*%U #should be equal to I
#dim(B%*%Bcoef)

g = B%*%Bcoef

#plot(1:nrow(g), g[,3], type="l")
#plot(g)
#matplot(g)
               
# Proportion of Variance Explained
pve = rho / sum(rho)
cum_pve = cumsum(pve)
               
# Plotting Functional Principal Components
n_plots <- ncol(g)
n_rows <- ceiling(n_plots / 2)   # 2 columns

#png("pca_bsplines.png", width = 1000, height = 400* n_rows)

par(mfrow = c(n_rows, 2))

for (i in 1:n_plots) {
  plot(t, g[, i], type = "l", lwd = 2, col = "blue",
       main = paste("Functional PC", i),
       xlab = "Time", ylab = "Value")
}

#dev.off()
               
# Reset graphics devices
while (!is.null(dev.list())) dev.off()

plot.new()
plot.window(xlim = c(0,1), ylim = c(0,1))
title("Test plot")

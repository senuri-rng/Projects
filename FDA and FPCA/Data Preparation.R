# Functional Data Analysis and FPCA
# Project 1: Canadian Precipitation Data

# Load data
data <- read.csv("Canadian_climate_history.csv")

# Select precipitation variables
precipitation_data <- data[, c(3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27)]

# Remove observations with missing values
Y <- na.omit(precipitation_data)

# Dimensions of the functional data
P <- nrow(Y)   # Number of time points
K <- ncol(Y)   # Number of weather stations

# Time domain
t <- seq(0, 1, length.out = P)

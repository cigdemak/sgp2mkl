require(dplyr) || install.packages("dplyr") 
library(dplyr)

source("../sgp/sgp_train.R")
source("../sgp/sgp_test.R")
source("../sgp2mkl/sgp2mkl_train.R")
source("../sgp2mkl/sgp2mkl_test.R")

pdist <- function(X1, X2) {
  if (identical(X1, X2) == TRUE) {
    D <- as.matrix(dist(X1))
  }
  else {
    D <- as.matrix(dist(rbind(X1, X2)))
    D <- D[1:nrow(X1), (nrow(X1) + 1):(nrow(X1) + nrow(X2))]
  }
  return(D)
}

#########################################################################################
###### DATA_NASA
#########################################################################################
data(nasa)
temperature <- matrix(0, 576, 72)
for (year in 1:6) {
  for (month in 1:12) {
    temperature[,(year - 1) * 12 + month] <- nasa$mets$surftemp[,,month,year]
  }
}
temperature <- temperature 

coordinates <- matrix(0, 576, 2)
colnames(coordinates) <- c("lat", "long")
coordinates[, "lat"] <- rep(nasa$dims$lat, 24)
coordinates[, "long"] <- rep(nasa$dims$long, each = 24)

periods <- matrix(0, 72, 2)
colnames(periods) <- c("year", "month")
periods[, "year"] <- rep(c(1995:2000), each = 12)
periods[, "month"] <- rep(1:12, 6)

cases <- temperature

X_temporal <- periods

#########################################################################################
###### DATA_NASA
#########################################################################################

s_ratio_latitude <- sqrt(0.5)
s_ratio_longitude <- sqrt(0.5)
s_ratio_year <- sqrt(0.5)
s_ratio_month <- sqrt(0.5)
parameters_fix <- c(s_ratio_latitude, s_ratio_longitude, s_ratio_year, s_ratio_month)

replication_number <- 1
sgp2mkl_spatiotemporal_results <- matrix(NA, replication_number, 2)
colnames(sgp2mkl_spatiotemporal_results) <- c("PCC", "NRMSE")
p <- matrix(NA, replication_number, 7)

for(repl in 1:replication_number){
  
  set.seed(1980 *repl)
  spatial_train_indices <- sample(1:nrow(cases), (nrow(cases)/2))
  spatial_test_indices <- setdiff(1:nrow(cases), spatial_train_indices)
  temporal_train_indices <- 1:floor(ncol(cases) * 5 / 6)
  temporal_test_indices <- floor(ncol(cases) * 5 / 6 + 1):ncol(cases)
  
  ## spatial data
  X_spatial_tra <- cbind(coordinates[spatial_train_indices, 1], coordinates[spatial_train_indices, 2])
  X_spatial_test <- cbind(coordinates[spatial_test_indices, 1], coordinates[spatial_test_indices, 2])
  
  X_spatial_tra <- scale(X_spatial_tra)
  X_spatial_test <- (X_spatial_test - matrix(attr(X_spatial_tra, "scaled:center"), nrow = nrow(X_spatial_test), ncol = ncol(X_spatial_test), byrow = TRUE)) / matrix(attr(X_spatial_tra, "scaled:scale"), nrow = nrow(X_spatial_test), ncol = ncol(X_spatial_test), byrow = TRUE)
  
  ## temporal data
  X_temporal_tra <- X_temporal[temporal_train_indices,]
  X_temporal_test <- X_temporal[temporal_test_indices,]
  
  D_longitude_tra_tra <- pdist(X_spatial_tra[,1, drop = FALSE], X_spatial_tra[,1, drop = FALSE])
  D_longitude_test_tra <- pdist(X_spatial_test[,1, drop = FALSE], X_spatial_tra[,1, drop = FALSE])
  D_longitude_test_test <- pdist(X_spatial_test[,1, drop = FALSE], X_spatial_test[,1, drop = FALSE])
  s <- s_ratio_longitude * mean(D_longitude_tra_tra)
  K_longitude_tra_tra <- exp(-D_longitude_tra_tra^2 / (2 * s^2))
  K_longitude_test_tra <- exp(-D_longitude_test_tra^2 / (2 * s^2))
  K_longitude_test_test <- exp(-D_longitude_test_test^2 / (2 * s^2))
  
  D_latitude_tra_tra <- pdist(X_spatial_tra[,2, drop = FALSE], X_spatial_tra[,2, drop = FALSE])
  D_latitude_test_tra <- pdist(X_spatial_test[,2, drop = FALSE], X_spatial_tra[,2, drop = FALSE])
  D_latitude_test_test <- pdist(X_spatial_test[,2, drop = FALSE], X_spatial_test[,2, drop = FALSE])
  s <- s_ratio_latitude * mean(D_latitude_tra_tra)
  K_latitude_tra_tra <- exp(-D_latitude_tra_tra^2 / (2 * s^2))
  K_latitude_test_tra <- exp(-D_latitude_test_tra^2 / (2 * s^2))
  K_latitude_test_test <- exp(-D_latitude_test_test^2 / (2 * s^2))
  
  D_year_tra_tra <- pdist(X_temporal_tra[,1, drop = FALSE], X_temporal_tra[,1, drop = FALSE])
  D_year_test_tra <- pdist(X_temporal_test[,1, drop = FALSE], X_temporal_tra[,1, drop = FALSE])
  D_year_test_test <- pdist(X_temporal_test[,1, drop = FALSE], X_temporal_test[,1, drop = FALSE])
  s <- s_ratio_year * mean(D_year_tra_tra)
  K_year_tra_tra <- exp(-D_year_tra_tra^2 / (2 * s^2))
  K_year_test_tra <- exp(-D_year_test_tra^2 / (2 * s^2))
  K_year_test_test <- exp(-D_year_test_test^2 / (2 * s^2))
  
  D_month_tra_tra <- pdist(X_temporal_tra[,2, drop = FALSE], X_temporal_tra[,2, drop = FALSE])
  D_month_test_tra <- pdist(X_temporal_test[,2, drop = FALSE], X_temporal_tra[,2, drop = FALSE])
  D_month_test_test <- pdist(X_temporal_test[,2, drop = FALSE], X_temporal_test[,2, drop = FALSE])
  s <- s_ratio_month * mean(D_month_tra_tra)
  K_month_tra_tra <- exp(-D_month_tra_tra^2 / (2 * s^2))
  K_month_test_tra <- exp(-D_month_test_tra^2 / (2 * s^2))
  K_month_test_test <- exp(-D_month_test_test^2 / (2 * s^2))
  
  Y_tra <- cases[spatial_train_indices, temporal_train_indices]
  mu <- mean(Y_tra)
  sig <- sd(Y_tra)
  Y_tra <- (Y_tra - mu)/sig
  
  Y_test <- cases[spatial_test_indices, temporal_test_indices]
  Y_test <- (Y_test - mu)/sig
  
  D_spatial_tra_tra <- array(c(D_latitude_tra_tra, D_longitude_tra_tra), c(nrow(D_latitude_tra_tra), ncol(D_latitude_tra_tra),2))
  D_temporal_tra_tra <- array(c(D_year_tra_tra, D_month_tra_tra), c(nrow(D_year_tra_tra), ncol(D_year_tra_tra),2))
  
  p[repl, ] <- sgp2mkl_train(D_temporal_tra_tra, D_spatial_tra_tra, Y_tra, parameters_fix)
  
  K_spatial_tra_tra <- p[2] * K_latitude_tra_tra + p[3] * K_longitude_tra_tra + p[4] * K_latitude_tra_tra * K_longitude_tra_tra
  K_spatial_test_tra <- p[2] * K_latitude_test_tra + p[3] * K_longitude_test_tra + p[4] * K_latitude_test_tra * K_longitude_test_tra
  K_spatial_test_test <- p[2] * K_latitude_test_test + p[3] * K_longitude_test_test + p[4] * K_latitude_test_test * K_longitude_test_test
  
  K_temporal_tra_tra <- p[5] * K_year_tra_tra + p[6] * K_month_tra_tra + p[7] * K_year_tra_tra * K_month_tra_tra
  K_temporal_test_tra <- p[5] * K_year_test_tra + p[6] * K_month_test_tra + p[7] * K_year_test_tra * K_month_test_tra
  K_temporal_test_test <- p[5] * K_year_test_test + p[6] * K_month_test_test + p[7] * K_year_test_test * K_month_test_test
  
  sigma <- p[repl, 1]
  parameters <- list(sigma = sigma)
  
  state <- sgp_train(K_spatial_tra_tra, K_temporal_tra_tra, Y_tra, parameters)
  prediction <- sgp2mkl_test(K_spatial_test_tra, K_spatial_test_test, K_temporal_test_tra, K_temporal_test_test, state)
  
  observed <- Y_test 
  predicted <- prediction$Y_mean 
  sgp2mkl_spatiotemporal_results[repl,] <- c(cor(matrix(observed, ncol = 1), matrix(predicted, ncol = 1)), sqrt(mean((observed - predicted)^2) / mean((observed - mean(observed))^2)))
}
print(sgp2mkl_spatiotemporal_results)
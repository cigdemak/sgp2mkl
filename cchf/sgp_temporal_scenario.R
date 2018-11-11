load("../data/coordinates.RData")
load("../data/case_counts.RData")

source("../sgp/sgp_train.R")
source("../sgp/sgp_test.R")

seasons <- c("1" = 1, "2" = 1, "3" = 1, "4" = 2, "5" = 3, "6" = 3, "7" = 3, "8" = 2, "9" = 2, "10" = 1, "11" = 1, "12" = 1)

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

s_ratio_latitude <- 0.25
s_ratio_longitude <- 0.25
s_ratio_year <- 1
s_ratio_month <- 1
s_ratio_season <- 1

set.seed(1980)
spatial_train_indices <- 1:nrow(cases)
spatial_test_indices <- 1:nrow(cases)
temporal_train_indices <- 1:floor(ncol(cases) * 10 / 12)
temporal_test_indices <- floor(ncol(cases) * 10 / 12 + 1):ncol(cases)

## spatial data
X_spatial_tra <- coordinates[spatial_train_indices, ]
rownames(X_spatial_tra) <- rownames(coordinates)[spatial_train_indices]
X_spatial_test <- coordinates[spatial_test_indices, ]
rownames(X_spatial_test) <- rownames(coordinates)[spatial_test_indices]

X_spatial_tra <- scale(X_spatial_tra)
X_spatial_test <- (X_spatial_test - matrix(attr(X_spatial_tra, "scaled:center"), nrow = nrow(X_spatial_test), ncol = ncol(X_spatial_test), byrow = TRUE)) / matrix(attr(X_spatial_tra, "scaled:scale"), nrow = nrow(X_spatial_test), ncol = ncol(X_spatial_test), byrow = TRUE)

## temporal data
X_temporal <- matrix(0, nrow = ncol(cases), ncol = 3)
for (month in 1:ncol(cases)) {
  X_temporal[month, 1] <- as.numeric(strsplit(colnames(cases)[month], "/")[[1]][1])
  X_temporal[month, 2] <- as.numeric(strsplit(colnames(cases)[month], "/")[[1]][2])
}
X_temporal[, 3] <- seasons[X_temporal[, 2]]

rownames(X_temporal) <- colnames(cases)

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
s <-  s_ratio_year * mean(D_year_tra_tra)
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

D_season_tra_tra <- pdist(X_temporal_tra[,3, drop = FALSE], X_temporal_tra[,3, drop = FALSE])
D_season_test_tra <- pdist(X_temporal_test[,3, drop = FALSE], X_temporal_tra[,3, drop = FALSE])
D_season_test_test <- pdist(X_temporal_test[,3, drop = FALSE], X_temporal_test[,3, drop = FALSE])
s <- s_ratio_season * mean(D_season_tra_tra)
K_season_tra_tra <- exp(-D_season_tra_tra^2 / (2 * s^2))
K_season_test_tra <- exp(-D_season_test_tra^2 / (2 * s^2))
K_season_test_test <- exp(-D_season_test_test^2 / (2 * s^2))

Y_tra <- log2(as.matrix(cases[spatial_train_indices, temporal_train_indices]) + 1)
Y_test <- log2(as.matrix(cases[spatial_test_indices, temporal_test_indices]) + 1)

p <- c(sd(Y_tra), rep(0, 2), rep(1,7))

K_spatial_tra_tra <- p[2] * K_latitude_tra_tra + p[3] * K_longitude_tra_tra + p[4] * K_latitude_tra_tra * K_longitude_tra_tra
K_spatial_test_tra <- p[2] * K_latitude_test_tra + p[3] * K_longitude_test_tra + p[4] * K_latitude_test_tra * K_longitude_test_tra
K_spatial_test_test <- p[2] * K_latitude_test_test + p[3] * K_longitude_test_test + p[4] * K_latitude_test_test * K_longitude_test_test

K_temporal_tra_tra <- p[5] * K_year_tra_tra + p[6] * K_month_tra_tra + p[7] * K_season_tra_tra + p[8] * K_year_tra_tra * K_month_tra_tra + p[9] * K_year_tra_tra * K_season_tra_tra + p[10] * K_month_tra_tra * K_season_tra_tra 
K_temporal_test_tra <- p[5] * K_year_test_tra + p[6] * K_month_test_tra + p[7] * K_season_test_tra + p[8] * K_year_test_tra * K_month_test_tra + p[9] * K_year_test_tra * K_season_test_tra + p[10] * K_month_test_tra * K_season_test_tra 
K_temporal_test_test <- p[5] * K_year_test_test + p[6] * K_month_test_test + p[7] * K_season_test_test + p[8] * K_year_test_test * K_month_test_test + p[9] * K_year_test_test * K_season_test_test + p[10] * K_month_test_test * K_season_test_test

sigma <- p[1]
parameters <- list(sigma = sigma)

state <- sgp_train(K_spatial_tra_tra, K_temporal_tra_tra, Y_tra, parameters)

prediction <- sgp_test(K_spatial_test_tra, K_spatial_test_test, K_temporal_test_tra, K_temporal_test_test, state)
observed <- 2^Y_test - 1
predicted <- 2^prediction$Y_mean - 1 
sgp_temporal_results <- c(cor(matrix(observed, ncol = 1), matrix(predicted, ncol = 1)), sqrt(mean((observed - predicted)^2) / mean((observed - mean(observed))^2)))
sprintf("PCC: %f, NRMSE: %f", sgp_temporal_results[1], sgp_temporal_results[2])
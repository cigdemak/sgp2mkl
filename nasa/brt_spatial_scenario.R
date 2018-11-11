total_replication_number <- 1

require(dismo) || install.packages("dismo") 
require(gbm) || install.packages("gbm") 
require(dplyr) || install.packages("dplyr") 
library(dismo)
library(gbm)
library(dplyr)

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

brt_spatial_results <- matrix(NA, total_replication_number, 2)
colnames(brt_spatial_results) <- c("PCC", "NRMSE")

for(repl in 1:total_replication_number){
  
  set.seed(1980 * repl)
  spatial_train_indices <- sample(1:nrow(cases), (nrow(cases)/2))
  spatial_test_indices <- setdiff(1:nrow(cases), spatial_train_indices)
  temporal_train_indices <- 1:ncol(cases) 
  temporal_test_indices <- 1:ncol(cases)  
  
  ## spatial data
  X_spatial_tra <- coordinates[spatial_train_indices,]
  X_spatial_test <- coordinates[spatial_test_indices,]
  X_spatial_tra <- scale(X_spatial_tra)
  X_spatial_test <- (X_spatial_test - matrix(attr(X_spatial_tra, "scaled:center"), nrow = nrow(X_spatial_test), ncol = ncol(X_spatial_test), byrow = TRUE)) / matrix(attr(X_spatial_tra, "scaled:scale"), nrow = nrow(X_spatial_test), ncol = ncol(X_spatial_test), byrow = TRUE)
  
  ## temporal data
  X_temporal_tra <- X_temporal[temporal_train_indices,]
  X_temporal_test <- X_temporal[temporal_test_indices,]
  
  data_tra <- as.data.frame(matrix(NA, nrow = nrow(X_spatial_tra) * nrow(X_temporal_tra), ncol = 5))
  colnames(data_tra) <- c("incidence", "year", "month", "latitude", "longitude")
  
  data_tra[,"incidence"] <- matrix((as.matrix(cases[spatial_train_indices, temporal_train_indices])), ncol = 1)
  data_tra[,"incidence"] <- scale(data_tra[,"incidence"])
  data_tra[,"year"] <- matrix(matrix(X_temporal_tra[,1], nrow = length(spatial_train_indices), ncol = length(temporal_train_indices), byrow = TRUE), ncol = 1)
  data_tra[,"month"] <- matrix(matrix(X_temporal_tra[,2], nrow = length(spatial_train_indices), ncol = length(temporal_train_indices), byrow = TRUE), ncol = 1)
  data_tra[,"latitude"] <- matrix(matrix(X_spatial_tra[,1], nrow = length(spatial_train_indices), ncol = length(temporal_train_indices), byrow = FALSE), ncol = 1)
  data_tra[,"longitude"] <- matrix(matrix(X_spatial_tra[,2], nrow = length(spatial_train_indices), ncol = length(temporal_train_indices), byrow = FALSE), ncol = 1)
  
  data_test <- as.data.frame(matrix(NA, nrow = nrow(X_spatial_test) * nrow(X_temporal_test), ncol = 5))
  colnames(data_test) <- c("incidence", "year", "month", "latitude", "longitude")
  
  data_test[,"incidence"] <- matrix((as.matrix(cases[spatial_test_indices, temporal_test_indices])), ncol = 1)
  data_test[,"incidence"] <- (data_test[,"incidence"] - matrix(attr(data_tra[,"incidence"], "scaled:center"), nrow = nrow(data_test[,"incidence"]), ncol = ncol(data_test[,"incidence"]), byrow = TRUE)) / matrix(attr(data_tra[,"incidence"], "scaled:scale"), nrow = nrow(data_test[,"incidence"]), ncol = ncol(data_test[,"incidence"]), byrow = TRUE)
  
  data_test[,"year"] <- matrix(matrix(X_temporal_test[,1], nrow = length(spatial_test_indices), ncol = length(temporal_test_indices), byrow = TRUE), ncol = 1)
  data_test[,"month"] <- matrix(matrix(X_temporal_test[,2], nrow = length(spatial_test_indices), ncol = length(temporal_test_indices), byrow = TRUE), ncol = 1)
  data_test[,"latitude"] <- matrix(matrix(X_spatial_test[,1], nrow = length(spatial_test_indices), ncol = length(temporal_test_indices), byrow = FALSE), ncol = 1)
  data_test[,"longitude"] <- matrix(matrix(X_spatial_test[,2], nrow = length(spatial_test_indices), ncol = length(temporal_test_indices), byrow = FALSE), ncol = 1)
  
  state <- gbm(incidence ~ year + month + latitude + longitude, data = data_tra, n.cores = 1,
               distribution = "gaussian", cv.folds = 5, n.trees = 100000, interaction.depth = 2)
  n_trees <- gbm.perf(state, method = "cv")
  prediction <- predict(state, data_test, n.trees = n_trees)
  
  predicted <- matrix(prediction, length(spatial_test_indices), length(temporal_test_indices))
  observed <- matrix(data_test$incidence, length(spatial_test_indices), length(temporal_test_indices))
  colnames(predicted) <- rownames(X_temporal_test)
  rownames(predicted) <- rownames(X_spatial_test)
  colnames(observed) <- rownames(X_temporal_test)
  rownames(observed) <- rownames(X_spatial_test)
  brt_spatial_results[repl,] <- c(cor(matrix(observed, ncol = 1), matrix(predicted, ncol = 1)), sqrt(mean((observed - predicted)^2) / mean((observed - mean(observed))^2)))
}
print(brt_spatial_results)
total_replication_number <- 1

require(randomForest) || install.packages("randomForest") 
library(randomForest)

load("../data/coordinates.RData")
load("../data/case_counts.RData")

seasons <- c("1" = 1, "2" = 1, "3" = 1, "4" = 2, "5" = 3, "6" = 3, "7" = 3, "8" = 2, "9" = 2, "10" = 1, "11" = 1, "12" = 1)

rf_spatial_results <-  matrix(NA, total_replication_number, 2)
colnames(rf_spatial_results) <- c("PCC", "NRMSE")

for(repl in 1:total_replication_number){
  
  set.seed(1980 * repl)
  spatial_train_indices <- sample(1:nrow(cases), 41)
  spatial_test_indices <- setdiff(1:nrow(cases), spatial_train_indices)
  temporal_train_indices <- 1:ncol(cases)
  temporal_test_indices <- 1:ncol(cases)
  
  ## spatial data
  X_spatial_tra <- coordinates[spatial_train_indices,]
  X_spatial_test <- coordinates[spatial_test_indices,]
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
  
  data_tra <- as.data.frame(matrix(NA, nrow = nrow(X_spatial_tra) * nrow(X_temporal_tra), ncol = 6))
  colnames(data_tra) <- c("incidence", "year", "month", "season", "latitude", "longitude")
  
  data_tra[,"incidence"] <- matrix(log2(as.matrix(cases[spatial_train_indices, temporal_train_indices]) + 1), ncol = 1)
  data_tra[,"year"] <- matrix(matrix(X_temporal_tra[,1], nrow = length(spatial_train_indices), ncol = length(temporal_train_indices), byrow = TRUE), ncol = 1)
  data_tra[,"month"] <- matrix(matrix(X_temporal_tra[,2], nrow = length(spatial_train_indices), ncol = length(temporal_train_indices), byrow = TRUE), ncol = 1)
  data_tra[,"season"] <- matrix(matrix(X_temporal_tra[,3], nrow = length(spatial_train_indices), ncol = length(temporal_train_indices), byrow = TRUE), ncol = 1)
  data_tra[,"latitude"] <- matrix(matrix(X_spatial_tra[,1], nrow = length(spatial_train_indices), ncol = length(temporal_train_indices), byrow = FALSE), ncol = 1)
  data_tra[,"longitude"] <- matrix(matrix(X_spatial_tra[,2], nrow = length(spatial_train_indices), ncol = length(temporal_train_indices), byrow = FALSE), ncol = 1)
  
  data_test <- as.data.frame(matrix(NA, nrow = nrow(X_spatial_test) * nrow(X_temporal_test), ncol = 6))
  colnames(data_test) <- c("incidence", "year", "month", "season", "latitude", "longitude")
  
  data_test[,"incidence"] <- matrix(log2(as.matrix(cases[spatial_test_indices, temporal_test_indices]) + 1), ncol = 1)
  data_test[,"year"] <- matrix(matrix(X_temporal_test[,1], nrow = length(spatial_test_indices), ncol = length(temporal_test_indices), byrow = TRUE), ncol = 1)
  data_test[,"month"] <- matrix(matrix(X_temporal_test[,2], nrow = length(spatial_test_indices), ncol = length(temporal_test_indices), byrow = TRUE), ncol = 1)
  data_test[,"season"] <- matrix(matrix(X_temporal_test[,3], nrow = length(spatial_test_indices), ncol = length(temporal_test_indices), byrow = TRUE), ncol = 1)
  data_test[,"latitude"] <- matrix(matrix(X_spatial_test[,1], nrow = length(spatial_test_indices), ncol = length(temporal_test_indices), byrow = FALSE), ncol = 1)
  data_test[,"longitude"] <- matrix(matrix(X_spatial_test[,2], nrow = length(spatial_test_indices), ncol = length(temporal_test_indices), byrow = FALSE), ncol = 1)
  
  state <- randomForest(incidence ~ year + month + season + latitude + longitude, data = data_tra, ntree = 100000)
  prediction <- predict(state, data_test)
  
  predicted <- matrix(2^prediction - 1, length(spatial_test_indices), length(temporal_test_indices))
  observed <- matrix(2^data_test$incidence - 1, length(spatial_test_indices), length(temporal_test_indices))
  colnames(predicted) <- rownames(X_temporal_test)
  rownames(predicted) <- rownames(X_spatial_test)
  colnames(observed) <- rownames(X_temporal_test)
  rownames(observed) <- rownames(X_spatial_test)
  rf_spatial_results[repl,] <- c(cor(matrix(observed, ncol = 1), matrix(predicted, ncol = 1)), sqrt(mean((observed - predicted)^2) / mean((observed - mean(observed))^2)))
}
print(rf_spatial_results)
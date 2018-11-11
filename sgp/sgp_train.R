sgp_train <- function(K_spatial_tra_tra, K_temporal_tra_tra, Y_tra, parameters) {
  
  svd_1 <- eigen(K_spatial_tra_tra)
  svd_2 <- eigen(K_temporal_tra_tra)
  
  state <- list(U_1 = svd_1$vectors, d_1 = svd_1$values, U_2 = svd_2$vectors, d_2 = svd_2$values, Y_tra = Y_tra, parameters = parameters)
}
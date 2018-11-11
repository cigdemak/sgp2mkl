sgp2mkl_train <- function(D_temporal_tra_tra, D_spatial_tra_tra, Y_tra, parameters_fix){
  n_spatial_kernels <- dim(D_spatial_tra_tra)[3] + choose(dim(D_spatial_tra_tra)[3], 2)
  n_temporal_kernels <- dim(D_temporal_tra_tra)[3] + choose(dim(D_temporal_tra_tra)[3], 2)
  
  gradient <- function(D_temporal_tra_tra, D_spatial_tra_tra, Y_tra, parameters_fix, parameters) {
    k_s <- array(NA, c(nrow(D_spatial_tra_tra), ncol(D_spatial_tra_tra), dim(D_spatial_tra_tra)[3]))
    k_t <- array(NA, c(nrow(D_temporal_tra_tra), ncol(D_temporal_tra_tra), dim(D_temporal_tra_tra)[3]))
    
    s_ratio_spatial <- as.numeric(parameters_fix[1:(dim(D_spatial_tra_tra)[3])])
    s_ratio_temporal <- as.numeric(tail(parameters_fix, dim(D_temporal_tra_tra)[3]))
    parameters_sigma <- as.numeric(parameters[1])
    
    interaction_spatial <- combn(dim(D_spatial_tra_tra)[3], 2)
    interaction_temporal <- combn(dim(D_temporal_tra_tra)[3], 2)
    D_spatial <- array(NA, ncol(interaction_spatial))
    D_temporal <- array(NA, (ncol(interaction_temporal)))
    
    k_s_interaction <- array(NA, c(nrow(D_spatial_tra_tra), ncol(D_spatial_tra_tra), ncol(interaction_spatial)))
    k_t_interaction <- array(NA, c(nrow(D_temporal_tra_tra), ncol(D_temporal_tra_tra), ncol(interaction_temporal)))
    
    parameters_spatial <- as.numeric(parameters[2:(1+(dim(D_spatial_tra_tra)[3] + ncol(interaction_spatial)))])
    parameters_temporal <- as.numeric(tail(parameters, (dim(D_temporal_tra_tra)[3] + ncol(interaction_temporal))))
    
    for(i in 1:dim(D_spatial_tra_tra)[3]) {
      k_s[,,i] <- exp((-D_spatial_tra_tra[ , ,i]^2) / (2 * (s_ratio_spatial[i] * mean(D_spatial_tra_tra[ , ,i]))^2))}
    for(i in 1:dim(D_temporal_tra_tra)[3]) {
      k_t[,,i] <- exp((-D_temporal_tra_tra[ , ,i]^2) / (2 * (s_ratio_temporal[i]* mean(D_temporal_tra_tra[ , ,i]))^2))}
    
    for(i in 1:ncol(interaction_spatial)) {
      k_s_interaction[,,i] <- apply(k_s[,,interaction_spatial[,i]], c(1,2), prod)}
    for(i in 1:ncol(interaction_temporal)) {
      k_t_interaction[,,i] <- apply(k_t[,,interaction_temporal[,i]], c(1,2), prod)}
    
    k_s <- abind(k_s, k_s_interaction)
    k_t <- abind(k_t, k_t_interaction)
    
    for(i in 1:length(parameters_spatial)) {
      k_s[,,i] <- parameters_spatial[i] * k_s[,,i]}
    for(i in 1:length(parameters_temporal)) {
      k_t[,,i] <- parameters_temporal[i] * k_t[,,i]}
    
    K_spatial_tra_tra <-  apply(k_s, c(1, 2), sum)
    K_temporal_tra_tra <- apply(k_t, c(1, 2), sum)
    
    N_1_tra <- ncol(K_spatial_tra_tra)
    N_2_tra <- ncol(K_temporal_tra_tra)
    
    state <- structured_gpr_train(K_spatial_tra_tra, K_temporal_tra_tra, Y_tra, list(sigma = parameters_sigma))
    
    diag_term <- (matrix(matrix(state$d_1, N_2_tra, N_1_tra, byrow = TRUE), N_1_tra * N_2_tra, 1) * matrix(matrix(state$d_2, N_2_tra, N_1_tra, byrow = FALSE), N_1_tra * N_2_tra, 1) + matrix(parameters_sigma^2, N_1_tra * N_2_tra, 1))^-1
    alpha <- matrix(state$U_2 %*% matrix(diag_term * (matrix(t(state$U_2) %*% t(state$Y_tra) %*% state$U_1, nrow = N_1_tra * N_2_tra, ncol = 1)), nrow = N_2_tra, ncol = N_1_tra) %*% t(state$U_1), nrow = N_1_tra * N_2_tra, ncol = 1) 
    D_sigma <-  (parameters_sigma * (t(alpha) %*% alpha)) - (1 / 2) * sum(2 * parameters_sigma * diag_term)
    
    for(s in 1:dim(k_s)[3]) {
      D_spatial[s] <- t(alpha) %*% (matrix(K_temporal_tra_tra %*% matrix(alpha, N_2_tra, N_1_tra) %*% t(k_s[,,s]), nrow = N_1_tra * N_2_tra, ncol = 1)) - 
        (1 / 2) * sum(diag_term * (matrix(matrix(diag(t(state$U_1) %*% k_s[,,s] %*% state$U_1), N_2_tra, N_1_tra, byrow = TRUE), N_1_tra * N_2_tra, 1) * 
                                     matrix(matrix(diag(t(state$U_2) %*% K_temporal_tra_tra %*% state$U_2), N_2_tra, N_1_tra, byrow = FALSE), N_1_tra * N_2_tra, 1)))
    }
    for (t in 1:dim(k_t)[3]) {
      D_temporal[t] <- t(alpha) %*% (matrix(k_t[,,t] %*% matrix(alpha, N_2_tra, N_1_tra) %*% t(K_spatial_tra_tra), nrow = N_1_tra * N_2_tra, ncol = 1)) - (1 / 2) * sum(diag_term * (matrix(matrix(diag(t(state$U_1) %*% K_spatial_tra_tra %*% state$U_1), N_2_tra, N_1_tra, byrow = TRUE), N_1_tra * N_2_tra, 1) * matrix(matrix(diag(t(state$U_2) %*% k_t[,,t] %*% state$U_2), N_2_tra, N_1_tra, byrow = FALSE), N_1_tra * N_2_tra, 1)))
    }
    
    return(c(-D_sigma, -D_spatial, -D_temporal))
  }
  
  loglikelihood <- function(D_temporal_tra_tra, D_spatial_tra_tra, Y_tra, parameters_fix, parameters) {
    k_t <- array(NA, c(nrow(D_temporal_tra_tra), ncol(D_temporal_tra_tra), dim(D_temporal_tra_tra)[3]))
    k_s <- array(NA, c(nrow(D_spatial_tra_tra), ncol(D_spatial_tra_tra), dim(D_spatial_tra_tra)[3]))
    
    s_ratio_spatial <- as.numeric(parameters_fix[1:(dim(D_spatial_tra_tra)[3])])
    s_ratio_temporal <- as.numeric(tail(parameters_fix, dim(D_temporal_tra_tra)[3]))
    parameters[parameters == 0] <- 1e-138
    parameters[parameters == Inf] <- 1e138
    parameters_sigma <- as.numeric(parameters[1])
    
    interaction_spatial <- combn(dim(D_spatial_tra_tra)[3], 2)
    interaction_temporal <- combn(dim(D_temporal_tra_tra)[3], 2)
    k_s_interaction <- array(NA, c(nrow(D_spatial_tra_tra), ncol(D_spatial_tra_tra), ncol(interaction_spatial)))
    k_t_interaction <- array(NA, c(nrow(D_temporal_tra_tra), ncol(D_temporal_tra_tra), ncol(interaction_temporal)))
    
    parameters_spatial <- as.numeric(parameters[2:(1+(dim(D_spatial_tra_tra)[3] + ncol(interaction_spatial)))])
    parameters_temporal <- as.numeric(tail(parameters, (dim(D_temporal_tra_tra)[3] + ncol(interaction_temporal))))
    
    for(i in 1:dim(D_spatial_tra_tra)[3]) {
      k_s[,,i] <- exp((-D_spatial_tra_tra[ , ,i]^2) / (2 * (s_ratio_spatial[i] * mean(D_spatial_tra_tra[ , ,i]))^2))}
    for(i in 1:dim(D_temporal_tra_tra)[3]) {
      k_t[,,i] <- exp((-D_temporal_tra_tra[ , ,i]^2) / (2 * (s_ratio_temporal[i]* mean(D_temporal_tra_tra[ , ,i]))^2))}
    
    for(i in 1:ncol(interaction_spatial)) {
      k_s_interaction[,,i] <- apply(k_s[,,interaction_spatial[,i]], c(1,2), prod)}
    for(i in 1:ncol(interaction_temporal)) {
      k_t_interaction[,,i] <- apply(k_t[,,interaction_temporal[,i]], c(1,2), prod)}
    
    k_s <- abind(k_s, k_s_interaction)
    k_t <- abind(k_t, k_t_interaction)
    
    for(i in 1:length(parameters_spatial)) {
      k_s[,,i] <- parameters_spatial[i] * k_s[,,i]}
    for(i in 1:length(parameters_temporal)) {
      k_t[,,i] <- parameters_temporal[i] * k_t[,,i]}
    
    K_spatial_tra_tra <- apply(k_s, c(1, 2), sum)
    K_temporal_tra_tra <- apply(k_t, c(1, 2), sum)
    
    N_1_tra <- ncol(K_spatial_tra_tra)
    N_2_tra <- ncol(K_temporal_tra_tra)
    
    state <- structured_gpr_train(K_spatial_tra_tra, K_temporal_tra_tra, Y_tra, list(sigma = parameters_sigma))
    diag_term <- (matrix(matrix(state$d_1, N_2_tra, N_1_tra, byrow = TRUE), N_1_tra * N_2_tra, 1) * matrix(matrix(state$d_2, N_2_tra, N_1_tra, byrow = FALSE), N_1_tra * N_2_tra, 1) + matrix(parameters_sigma^2, N_1_tra * N_2_tra, 1))^-1
    alpha <- matrix(state$U_2 %*% matrix(diag_term * (matrix(t(state$U_2) %*% t(state$Y_tra) %*% state$U_1, nrow = N_1_tra * N_2_tra, ncol = 1)), nrow = N_2_tra, ncol = N_1_tra) %*% t(state$U_1), nrow = N_1_tra * N_2_tra, ncol = 1) 
    diag_term_det <- (matrix(matrix(state$d_1, N_2_tra, N_1_tra, byrow = TRUE), N_1_tra * N_2_tra, 1) * matrix(matrix(state$d_2, N_2_tra, N_1_tra, byrow = FALSE), N_1_tra * N_2_tra, 1))
    diag_term_det[diag_term_det == Inf] <- 1e138
    
    sigma_sq <- parameters_sigma^-2
    sigma_sq[sigma_sq == Inf] <-1e138
    log_det <- sum(log((diag_term_det^-1 + sigma_sq) * (diag_term_det) * parameters_sigma^2))
    
    loglikelihood <- -(1 / 2) * (t(matrix(t(Y_tra), N_1_tra * N_2_tra, 1)) %*% alpha) - (1 / 2) * log_det - ((N_1_tra * N_2_tra) / 2) * log(2 * pi)
    loglikelihood[loglikelihood == -Inf] <- -1e138
    
    return(-loglikelihood)  
  }
  
  l <- function(x) {
    loglikelihood(D_temporal_tra_tra, D_spatial_tra_tra, Y_tra, parameters_fix, exp(x))
  }
  g <- function(x) {
    gradient(D_temporal_tra_tra, D_spatial_tra_tra, Y_tra, parameters_fix, exp(x))
  }
  
  heq <- function(x){
    ui <- matrix(c(0, rep(1, n_spatial_kernels), rep(0, n_temporal_kernels),
                   0, rep(0, n_spatial_kernels), rep(1, n_temporal_kernels)), 2, (n_spatial_kernels  + n_temporal_kernels + 1), byrow = TRUE)
    return(rowSums(ui * matrix(exp(x), 2, length(x), byrow = TRUE)) - c(1, 1))
  }
  
  heq.jac <- function(x){
    ui <- matrix(c(0, rep(1, n_spatial_kernels), rep(0, n_temporal_kernels),
                   0, rep(0, n_spatial_kernels), rep(1, n_temporal_kernels)), 2, (n_spatial_kernels + n_temporal_kernels + 1), byrow = TRUE)
    return(ui * matrix(exp(x), 2, length(x), byrow = TRUE))
  }
  
  res <- constrOptim.nl(par = log(c(1, rep(1/n_spatial_kernels, n_spatial_kernels), rep(1/n_temporal_kernels, n_temporal_kernels))), 
                        fn = l, 
                        gr = g, 
                        heq = heq, 
                        heq.jac = heq.jac)
  
  par_rep_temporal <- exp(res$par)
}
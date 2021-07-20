sample_Y_given_Z = function(Z, SNR, sd_eps){
  n = nrow(Z)
  p = ncol(Z)
  sd_gamma = sqrt(SNR*sd_eps^2/mean(rowSums(Z*Z)))
  gamma = matrix(rnorm(p, mean = 0, sd = sd_gamma),p,1)
  epsilon = matrix(rnorm(n, mean = 0, sd = sd_eps),n,1)
  Y = Z%*%gamma + epsilon
}

sample_X = function(n, pi_init){
  X = matrix(rbinom(n = n, size = 1, prob = pi_init), n, 1)
  return(X)
}

sample_Z = function(n, p, pi_init, pi_flip){
  X = sample_X(n, pi_init)
  Z = sample_Z_given_X(X, p, pi_flip)
  return(Z)
}

sample_Z_given_X = function(X, p, pi_flip){
  n = nrow(X)
  flip  = matrix(rbinom(n = n*p, size = 1, prob = pi_flip), n, p)
  Z = apply(t(apply(flip, 1, cumsum)), 2, function(col)(X+col)) %% 2
  return(Z)
}

mean_X_given_Z = function(Z, pi_init, pi_flip){
  prob_0 = (1-pi_init)*(pi_flip*Z[,1] + (1-pi_flip)*(1-Z[,1]))
  prob_1 = pi_init*(pi_flip*(1-Z[,1]) + (1-pi_flip)*Z[,1])
  prob_vec = prob_1/(prob_0 + prob_1) 
  return(prob_vec)
}

sample_X_given_Z = function(Z, pi_init, pi_flip){
  n = nrow(Z)
  prob_vec = mean_X_given_Z(Z, pi_init, pi_flip)
  X = matrix(rbinom(n = n, size = 1, prob = prob_vec),n,1)
  return(X)
}

var_X_given_Z = function(Z, pi_init, pi_flip){
  prob_vec = mean_X_given_Z(Z, pi_init, pi_flip)
  return(prob_vec*(1-prob_vec))
}

compute_U = function(X_test, Y_test, Z_test, g_hat, pi_init, pi_flip){
  n_test = nrow(X_test)
  g_hat_test = predict(g_hat, newx = Z_test)
  statistic = sum((Y_test - g_hat_test)*(X_test - mean_X_given_Z(Z_test, pi_init, pi_flip)))
  normalization = mean((Y_test - g_hat_test)^2*var_X_given_Z(Z_test, pi_init, pi_flip))
  U = statistic/sqrt(n_test*normalization)
  return(U)
}
library(MASS)
library(evd)
library(dplyr)
data = wind_speed
B = 24                # Block size
num_thetas = B/2-1    # Number of theta values
T = nrow(data) # Sample size
p = ncol(data) # Dimensions of time series

# Set of frequencies theta
Theta = 2*(1:num_thetas)*pi/B  # Set of frequencies

# Spectral density estimator
spectral_density_estimator = function(theta, X, B){
  # X: Time series data as an (T x p) matrix
  # B: Block size
  
  T = nrow(X)  # Sample size
  p = ncol(X)  # Dimension
  
  # Check that T is divisible by B
  if (T %% B != 0) {
    stop("The sample size T must be divisible by the block size B.")
  }
  
  # Initialize F_hat
  F_hat = matrix(0, ncol = p, nrow = p)  # p x p matrix
  
  # Loop over each block r
  for (r in 1:(T/B)) {
    
    # Indices for the current block
    start_idx = (r - 1)*B + 1
    end_idx = r*B
    
    # Extract the block of data
    X_block = X[start_idx:end_idx, ]
    
    # Compute the first summation
    sum_1 = colSums(X_block*exp(-1i*(start_idx:end_idx)*theta))
    
    # Compute the second summation
    sum_2 = colSums(X_block*exp(1i*(start_idx:end_idx)*theta))
    
    # Multiply the results and add to the final sum
    F_hat = F_hat + sum_1 %*% t(sum_2)
  }
  
  # Final estimator
  F_hat = F_hat/(2*pi*T)
  
  return(F_hat)
}


# Covariance matrix for each theta 
cov_theta = function(F_theta){
  p = 1
  
  # Initialize a covariance matrix Sigma
  Sigma = matrix(0, nrow = 2*p^2, ncol = 2*p^2)
  
  # Var(Re_12)
  Sigma[1,1] = 1/2 * Re(F_theta[1,1] * F_theta[2,2]) + 1/2 * Re(F_theta[1,2]^2)
  
  # Var(Im_12)
  Sigma[2,2] = 1/2 * Re(F_theta[1,1] * F_theta[2,2]) - 1/2 * Re(F_theta[1,2]^2)
  
  # Cov(Re_12, Im_12) and Cov(Im_12, Re_12)
  Sigma[1,2] = Sigma[2,1] = 1/2 * Im(F_theta[1,2]^2)
  
  return(Sigma)
}

# Function to compute test statistics and 95 quantile
test_statistic_and_quantile = function(data_pairs){
  
  # Initialize F_theta_hat and F_star
  F_theta_hat_list = list()
  F_star = c()
  
  # Compute F_theta_hat and F_star
  for (k in 1:num_thetas){
    F_theta_hat_list[[k]] = spectral_density_estimator(Theta[k], data_pairs, B)
    F_star[k] = max(abs(Re(F_theta_hat_list[[k]][1,2])), abs(Im(F_theta_hat_list[[k]][1,2])))
  }
  
  # Compute test_statistic
  test_statistic = max(sqrt(T/B) * F_star)
  
  # Initialize the theoretical covariance matrix Sigma_theoretical
  p = 1
  Sigma_theoretical = matrix(0, nrow = 2*num_thetas*p^2, ncol = 2*num_thetas*p^2)
  
  # Loop through each covariance matrix and place it in the diagonal
  for (i in 1:num_thetas){
    # Calculate the start and end indices for the diagonal block
    start_idx = (i - 1) * 2*p^2 + 1
    end_idx = i * 2*p^2
    
    
    # Place the covariance matrix in the diagonal block
    Sigma_theoretical[start_idx:end_idx, start_idx:end_idx] = block_list[[i]] =  cov_theta(F_theta_hat_list[[i]])
  }
  
  
  # Function to compute the max norm (infinity norm)
  max_norm = function(x){
    return(max(abs(x)))
  }
  
  # epsilon = 1e-6
  # Sigma_reg = Sigma_theoretical + diag(epsilon, nrow(Sigma))
  
  # Sigma_nearPD = nearPD(Sigma_theoretical)$mat
  
  # Simulate the distribution of |Z|_infty
  n_sim = 10000  # Number of simulations
  # Z_sim = mvrnorm(n_sim, mu = rep(0, (B-2)*p^2), Sigma = as.matrix(Sigma_nearPD))  # Simulate from N(0, Sigma)
  Z_sim = mvrnorm(n_sim, mu = rep(0, (B-2)*p^2), Sigma = Sigma_theoretical)  # Simulate from N(0, Sigma)
  max_values = apply(Z_sim, 1, max_norm)  # Compute the infinity norm for each row
  
  # Find the 95% cutoff value
  cutoff_95 = quantile(max_values, 0.95)
  
  return(c(test_statistic, cutoff_95))
}

results = list()

# Loop over all unique pairs of time series
for (a in 1:(p - 1)) {
  for (b in (a + 1):p) {
    
    data_pairs = data[, c(a, b)]
    stats = test_statistic_and_quantile(data_pairs)
    pair_name = paste0("Pair (", a, ", ", b, ")")
    
    results[[pair_name]] = data.frame(Pair = pair_name,
                                       Test_Statistic = stats[1],
                                       Quantile_95 = stats[2])
  }
}

# Combine all results into data frame
results_df = bind_rows(results)

# view the results
print(results_df)

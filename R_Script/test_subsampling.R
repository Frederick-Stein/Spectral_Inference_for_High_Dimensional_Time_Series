library(MASS)
library(evd)
library(dplyr)
data_wind = wind_speed
B = 24                 # Block size
num_thetas = B/2-1    # Number of theta values
T = nrow(data_wind) # Sample size
M = 288                # Subsample size
B_sub = 24            # Subsample Block size
num_thetas_sub = B_sub/2-1 # Number of theta values in subsample
p = ncol(data_wind) # Dimensions of time series

# Set of frequencies theta
Theta = 2*(1:num_thetas)*pi/B  # Set of frequencies
Theta_sub = 2*(1:num_thetas_sub)*pi/B_sub  # Set of frequencies



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





# Function to compute F_l
F_l = function(data_pairs, F_theta_hat_list, M, B_sub){
  
  T = nrow(data_pairs) # Sample size
  L = T - M + 1       # Number of F_l
  
  # Initialize vector of F_l
  F_l_list = c()
  
  for (l in 1:L){
    # Extract subsample
    data_pairs_sub =  data_pairs[l:(l+M-1),]
    
    
    f_l = 0
    
    # Compute F_l
    for (k in 1:num_thetas_sub){
      
      # Compute F_l_theta_hat
      F_l_theta_hat = spectral_density_estimator(Theta_sub[k], data_pairs_sub, B_sub)
      
      error = F_l_theta_hat[1,2] - F_theta_hat_list[[k]][1,2]
      
      max_Re = max(abs(Re(error)))
      max_Im = max(abs(Im(error)))
      
      F_l_theta = max(max_Re, max_Im)
      f_l = max(f_l,  F_l_theta)
    }
    
    # Store the l-th sample result in F_l_list
    F_l_list[l] = sqrt(M/B_sub) * f_l
  }
  return(F_l_list)
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
  
  # Compute F_l
  F_l_list = F_l(data_pairs, F_theta_hat_list, M, B_sub)
  
  # Find the 95% cutoff value
  cutoff_95 = quantile(F_l_list, 0.95)
  return(c(test_statistic, cutoff_95))
}

results = list()

# Loop over all unique pairs of time series
for (a in 1:(p - 1)) {
  for (b in (a + 1):p) {
    
    data_pairs = data_wind[, c(a, b)]
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

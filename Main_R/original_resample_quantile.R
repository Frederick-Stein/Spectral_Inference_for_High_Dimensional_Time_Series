library(MASS)
library(evd)
library(Matrix)

# Parameters
p = 15                # Dimension
T = 2160              # Sample size
rho = 0.5            # Correlation coefficient
N = 2001              # Number of observations
B = 30              # Block size
# M = 200              # Subsample size
# B_sub = 20           # Subsample block size
num_thetas = B/2-1   # Number of theta values
# num_thetas_sub = B_sub/2-1 # Number of theta values in subsample
# ratio = B/B_sub
alpha = 2            # Shape parameter of the Gamma distribution
beta = 1             # Rate parameter of the Gamma distribution
n_sim = 1000           # Number of simulations

Cov_true = alpha / beta^2 * diag(p)  # Theoretical covariance matrix

# Toeplitz matrix
tp = c(rho^(1:p))
Tp = toeplitz(tp)

# Create the list of matrices A_k for each lag k
A_list = lapply(0:(N-1), function(k) rho^k * Tp)

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

# Function to compute the covariance matrix for each theta
cov_theta = function(F_theta){
  p = nrow(F_theta)
  Sigma = matrix(0, nrow = 2*p^2, ncol = 2*p^2)
  
  for (i in 1:p) {
    for (j in 1:p) {
      for (h in 1:p) {
        for (l in 1:p) {
          # Covariances
          cov_Reij_Rehl = 0.5 * (Re(F_theta[i,h] * F_theta[l,j]) + Re(F_theta[i,l] * F_theta[h,j]))
          cov_Imij_Imhl = 0.5 * (Re(F_theta[i,h] * F_theta[l,j]) - Re(F_theta[i,l] * F_theta[h,j]))
          cov_Reij_Imhl = 0.5 * (Im(F_theta[i,l] * F_theta[h,j]) - Im(F_theta[i,h] * F_theta[l,j]))
          cov_Imij_Rehl = 0.5 * (Im(F_theta[i,l] * F_theta[h,j]) + Im(F_theta[i,h] * F_theta[l,j]))
          
          # Indices for filling Sigma
          idx_A1 = (i-1)*p+j
          idx_A2 = (h-1)*p+l
          idx_B1 = idx_A1 + p^2
          idx_B2 = idx_A2 + p^2
          
          # Fill Sigma matrix
          Sigma[idx_A1, idx_A2] = cov_Reij_Rehl
          Sigma[idx_B1, idx_B2] = cov_Imij_Imhl
          Sigma[idx_A1, idx_B2] = cov_Reij_Imhl
          Sigma[idx_B1, idx_A2] = cov_Imij_Rehl
        }
      }
    }
  }
  
  return(Sigma)
}

# Function to generate the process X
generate_process_X = function(T, N, p, Tp, A_list, alpha, beta){
  eps = matrix(rgamma((N+T-1)*p, shape = alpha, rate = beta), nrow = N+T-1, ncol = p) - alpha / beta
  X = matrix(0, nrow = T, ncol = p)
  
  for (k in 0:(N-1)) {
    X[1, ] = X[1, ] + A_list[[k+1]] %*% eps[N-k, ]
  }
  for (t in 1:(T-1)) {
    X[t+1, ] =  Tp %*% eps[t+N, ] + rho * (X[t, ] - A_list[[N]] %*% eps[t, ])
  }
  
  return(X)
}

# Function to compute theoretical spectral density
spectral_density_theoretical = function(theta, Tp, Cov_true, N){
  A_theta = matrix(0, nrow = nrow(Tp), ncol = ncol(Tp))
  
  # Fourier transform over lags
  for (k in 0:(N-1)) {
    A_theta = A_theta + A_list[[k+1]] * exp(-1i * theta * k)
  }
  
  # Spectral density
  f_theta = (1/(2*pi)) * A_theta %*% Cov_true %*% Conj(t(A_theta))
  return(f_theta)
}

# Function to generate Z_hat_values blockwise
generate_Z_hat_blockwise = function(Sigma_hat, num_samples) {
  block_sizes = sapply(Sigma_hat, nrow)  # Size of each block
  Z_hat_values = matrix(0, nrow = num_samples, ncol = sum(block_sizes))  # Initialize Z_hat_values
  
  # Generate samples for each block
  start_idx = 1
  for (b in seq_along(Sigma_hat)) {
    block_size = block_sizes[b]
    Z_hat_block = mvrnorm(num_samples, mu = rep(0, block_size), Sigma = Sigma_hat[[b]])  # Generate block samples
    
    # Fill the corresponding part of Z_hat_values
    Z_hat_values[, start_idx:(start_idx + block_size - 1)] = Z_hat_block
    start_idx = start_idx + block_size
  }
  
  return(Z_hat_values)
}

# Main simulation loop with blockwise Z_hat generation
run_simulation_blockwise = function(n_sim){
  Theta = 2 * (1:num_thetas) * pi / B
  Theta_sub = 2 * (1:num_thetas_sub) * pi / B_sub
  
  Z_hat_quantiles_95 = c()
  Z_hat_quantiles_90 = c()
  F_star_values = c()
  
  # Compute quantiles of |Z_hat|_infty
  for (n in 1:n_sim) {
    # Generate process X
    X = generate_process_X(T, N, p, Tp, A_list, alpha, beta)
    
    
    # Compute F_theta_hat_list
    F_theta_hat_list = lapply(1:num_thetas, function(i) {
      spectral_density_estimator(Theta[i], X, B)
    })
    
    # Covariance matrices and Z_hat_quantiles
    Sigma_hat_blocks = lapply(F_theta_hat_list, function(F_theta) {
      cov_theta(F_theta)
    })
    
    # Generate Z_hat values blockwise
    Z_hat_values = generate_Z_hat_blockwise(Sigma_hat_blocks, 1000)
    Z_hat_quantile_95 = quantile(apply(Z_hat_values, 1, function(x) max(abs(x))), 0.95)
    Z_hat_quantile_90 = quantile(apply(Z_hat_values, 1, function(x) max(abs(x))), 0.9)
    
    # Store Z_hat quantiles
    Z_hat_quantiles_95[n] = Z_hat_quantile_95
    Z_hat_quantiles_90[n] = Z_hat_quantile_90
    
    # Compute F_star
    F_star_max_real = 0
    F_star_max_imag = 0
    
    for (i in 1:num_thetas){
      F_theta = spectral_density_theoretical(Theta[i], Tp, Cov_true, N)
      error_real = max(abs(Re(F_theta_hat_list[[i]] - F_theta)))
      error_imag = max(abs(Im(F_theta_hat_list[[i]] - F_theta)))
      F_star_max_real = max(F_star_max_real, error_real)
      F_star_max_imag = max(F_star_max_imag, error_imag)
    }
    
    F_star = sqrt(T/B) * max(F_star_max_real, F_star_max_imag)
    F_star_values[n] = F_star
  }
  
  return(list(Z_hat_quantiles_95 = Z_hat_quantiles_95, Z_hat_quantiles_90 = Z_hat_quantiles_90, 
              F_star_values = F_star_values))
}

# Run the simulation with blockwise generation
simulation_results_blockwise = run_simulation_blockwise(n_sim)

F_star_Z_hat_95_blockwise = sum(simulation_results_blockwise$F_star_values < simulation_results_blockwise$Z_hat_quantiles_95)
F_star_Z_hat_90_blockwise = sum(simulation_results_blockwise$F_star_values < simulation_results_blockwise$Z_hat_quantiles_90)

print(paste("p =", p, "T =", T, "B =", B))
print(paste("F_star vs Z_hat 95% quantile: " , F_star_Z_hat_95_blockwise))
print(paste("F_star vs Z_hat 90% quantile: " , F_star_Z_hat_90_blockwise))

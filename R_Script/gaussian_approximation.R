library(MASS)
library(evd)
# Parameters
p = 15                # Dimension
T = 2160               # Sample size
rho = 0.5             # Correlation coefficient
N = 2001              # Number of observations
B = 24                 # Block size
num_thetas = B/2-1    # Number of theta values

# Cov_true = diag(p)  # Theoretical covariance matrix of normal-distributed epsilon

alpha = 2           # Shape parameter of the Gamma distribution
beta = 1              # Rate parameter of the Gamma distribution
Cov_true = alpha/beta^2 * diag(p)  # Theoretical covariance matrix of gamma-distributed epsilon

# loc = 0             # Location parameter (mu)
# scale = 2            # Scale parameter (beta)
# Cov_true = (pi^2/6) * scale^2 * diag(p)  # Theoretical covariance matrix of gumbel-distributed epsilon

# Number of simulations
n_sim = 1000

# Set of frequencies theta
Theta = 2*(1:num_thetas)*pi/B  # Set of frequencies

# Toeplitz matrix
tp = c(rho^(1:p))
Tp = toeplitz(tp)
# A = t(matrix(rep(c(Tp), N), p*p, N))*rho^(0:(N-1))
# A_matrices = array(t(A), dim = c(p, p, N))  # Reshape A into a 3D array (p x p x N)
# Create the list of matrices A_k for each lag k
# A_list = list()
# for (k in 0:(N-1)){
#   A_list[[k+1]] = rho^k * Tp
# }

A_list = list()
for (k in 0:(N-1)){
  A_list[[k+1]] = rho^k * Tp
}


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
  # Initialize a covariance matrix Sigma
  p = nrow(F_theta)
  Sigma = matrix(0, nrow = 2*p^2, ncol = 2*p^2)
  
  # # Extract the real and imaginary parts
  # Re_part = Re(F_hat_theta)
  # Im_part = Im(F_hat_theta)
  
  # Helper function to compute linear index from (i, j) pair for vectorized matrices
  index = function(i, j, p){
    return((i-1)*p+j)
  }
  
  # Compute covariance matrix
  for (i in 1:p){
    for (j in 1:p){
      for (h in 1:p){
        for (l in 1:p){
          
          # Cov(Re_ij, Re_hl)
          cov_Reij_Rehl = 0.5*(Re(F_theta[i, h]*F_theta[l, j]) + Re(F_theta[i, l]*F_theta[h, j]))
          
          # Cov(Im_ij, Im_hl)
          cov_Imij_Imhl = 0.5*(Re(F_theta[i, h]*F_theta[l, j]) - Re(F_theta[i, l]*F_theta[h, j]))
          
          # Cov(Re_ij, Im_hl)
          cov_Reij_Imhl = 0.5*(Im(F_theta[i, l]*F_theta[h, j]) - Im(F_theta[i, h]*F_theta[l, j]))
          
          # Cov(Im_ij, Re_hl)
          cov_Imij_Rehl = 0.5*(Im(F_theta[i, l]*F_theta[h, j]) + Im(F_theta[i, h]*F_theta[l, j]))
          
          # Get the linear indices
          idx_A1 = index(i, j, p)            
          idx_A2 = index(h, l, p)           
          idx_B1 = index(i, j, p) + p^2      
          idx_B2 = index(h, l, p) + p^2     
          
          # Fill the Sigma matrix with the computed covariances
          Sigma[idx_A1, idx_A2] = cov_Reij_Rehl      # Cov(Re_ij, Re_hl)
          Sigma[idx_B1, idx_B2] = cov_Imij_Imhl      # Cov(Im_ij, Im_hl)
          Sigma[idx_A1, idx_B2] = cov_Reij_Imhl      # Cov(Re_ij, Im_hl)
          Sigma[idx_B1, idx_A2] = cov_Imij_Rehl      # Cov(Im_ij, Re_hl)
        }
      }
    }
  }
  
  return(Sigma)
}

# Function to compute the theoretical spectral density
spectral_density_theoretical = function(theta, Tp, Cov_true, N){
  # Initialize the Fourier transform sum
  A_theta = matrix(0, nrow = nrow(Tp), ncol = ncol(Tp))
  
  # Loop over lags to compute the Fourier transform of A_k
  for (k in 0:(N-1)){
    A_theta = A_theta + A_list[[k+1]] * exp(-1i * theta * k)
  }
  
  # Compute the spectral density
  f_theta = (1/(2*pi)) * A_theta %*% Cov_true %*% Conj(t(A_theta))
  
  return(f_theta)
}

# Initialize the theoretical covariance matrix Sigma_theoretical
Sigma_theoretical = matrix(0, nrow = 2*num_thetas*p^2, ncol = 2*num_thetas*p^2)
# Initialize F_theta list
F_theta_list = list()
# Initialize the block covariance
block_list = list()


# Loop through each covariance matrix and place it in the diagonal
for (i in 1:num_thetas){
  # Calculate the start and end indices for the diagonal block
  start_idx = (i - 1) * 2*p^2 + 1
  end_idx = i * 2*p^2
  
  # Compute F_theta and place it in F_theta_list
  F_theta_list[[i]] = spectral_density_theoretical(Theta[i], Tp, Cov_true, N)
  
  # Place the covariance matrix in the diagonal block
  Sigma_theoretical[start_idx:end_idx, start_idx:end_idx] = block_list[[i]] = cov_theta(F_theta_list[[i]])
}





# Store F_star values
F_star_values = numeric(n_sim)

for (i in 1:n_sim){
  # Generate data X
  # Generate Gamma-distributed noise
  # eps = matrix(rgamma(N * p, shape = alpha, rate = beta), N, p)
  # eps = matrix(rnorm(N * p, 0, 1), N, p)
  # 
  # # Fourier transform operations
  # X = matrix(0, T, p)
  # z = mvfft(A) * matrix(Conj(mvfft(eps)), N, p * p)
  # x = Re(mvfft(z, inverse = TRUE))[(N - T + 1):N, ] / N
  # y = array(x, c(T, p, p))
  # 
  # # Summing over dimensions to get X
  # X = apply(y, c(1, 3), sum)
  
  
  
  # Generate Gamma-distributed noise epsilon
  eps = matrix(rgamma((N+T-1)*p, shape = alpha, rate = beta), nrow = N+T-1, ncol = p) - alpha/beta
  
  # Generate normal-distributed epsilon
  # eps = matrix(rnorm((N+T-1) * p, 0, 1), N+T-1, p)
  
  # Generate Gumbel-distributed epsilon
  # eps = matrix(rgumbel((N + T - 1) * p, loc = loc, scale = scale), nrow = N + T - 1, ncol = p) - (loc + scale * 0.5772)
  
  # Initialize X matrix to store the generated process
  X = matrix(0, nrow = T, ncol = p)  # X will store the generated data
  
  # # Compute the linear process X_i = sum(A_k * epsilon_{i-k})
  # for (t in 1:T) {
  #   for (k in 0:(N-1)){
  #     X[t, ] = X[t, ] + A_list[[k+1]] %*% eps[t+N-1-k, ]
  #   }
  # }
  
  # Compute the linear process X_i = sum(A_k * epsilon_{i-k})
  for (k in 0:(N-1)){
    X[1, ] = X[1, ] + A_list[[k+1]] %*% eps[N-k, ]
  }
  for (t in 1:(T-1)){
    X[t+1, ] = Tp %*% eps[t+N, ] + rho * (X[t, ] - A_list[[N]] %*% eps[t, ])
  }
  
  
  
  # Compute the empirical spectral density and F_star
  max_error_real = 0
  max_error_imag = 0
  
  for (k in 1:num_thetas){
    # Compute empirical spectral density for the current theta
    F_empirical = spectral_density_estimator(Theta[k], X, B)
    
    # Compute theoretical spectral density for the current theta
    F_theoretical = F_theta_list[[k]]
    
    # Compute the difference
    error = F_empirical - F_theoretical
    
    # Compute the infinity norms of the real and imaginary parts
    error_real_norm = max(abs(Re(error)))
    error_imag_norm = max(abs(Im(error)))
    
    # Track the maximum norms across all theta
    max_error_real = max(max_error_real, error_real_norm)
    max_error_imag = max(max_error_imag, error_imag_norm)
  }
  
  # Compute F_star
  F_star_values[i] = sqrt(T/B) * max(max_error_real, max_error_imag)
}



# Generate |Z|_infty

# Function to compute the max norm (infinity norm)
max_norm = function(x){
  return(max(abs(x)))
}


# Generate data block by block and combine
block_samples_list = list()

for (i in 1:num_thetas) {
  block_cov = block_list[[i]]
  
  # Generate multivariate normal data for this block
  block_samples = mvrnorm(n_sim, mu = rep(0, ncol(block_cov)), Sigma = block_cov)
  
  # Store the samples in a list
  block_samples_list[[i]] = block_samples
}

# Combine the samples from all blocks
results_true = apply(do.call(cbind, block_samples_list), 1, max_norm)

# # Generate theoretical Z with covariance Sigma_theoretical
# Mean = rep(0, (B-2)*p^2)
# Z = mvrnorm(n_sim, mu = Mean, Sigma = Sigma_theoretical)
# 
# results_true = numeric(n_sim)
# 
# for (i in 1:n_sim){
#   # Compute the infinity norm
#   results_true[i] = max_norm(Z[i,])
# }


# QQ-plot comparing the simulated data
qqplot(results_true, F_star_values, 
       # main = "QQ-Plot: Estimated vs True",
       xlab = "Gaussian Approximation", ylab = "Spectral Estimator", cex = 0.8)
abline(0, 1, col = "red")

# Add a text annotation in the top left corner
text(x = min(par("usr")[1:2]),  # X coordinate: min X value of the plot
     y = max(par("usr")[3:4]) - 0.1 * diff(par("usr")[3:4]),  # Lower Y coordinate
     labels = paste("p =", p, " T =", T, "\nB =", B, " n =", n_sim),  # Annotation text
     pos = 4,  # Position: 4 = to the right of (x, y)
     cex = 1,  # Text size
     col = "blue",  # Text color
     adj = c(0, 1))  # Fine-tuning text position


# Calculate quantiles
results_true_quantiles = quantile(results_true, probs = ppoints(length(results_true)))
F_star_values_quantiles = quantile(F_star_values, probs = ppoints(length(F_star_values)))

# Fit a linear model to the quantiles
fit = lm(F_star_values_quantiles ~ results_true_quantiles)
print(fit)

# Load necessary libraries
# install.packages("MASS") # If not already installed
library(MASS) # For mvrnorm

# Helper function to calculate Euclidean distance matrix
calculate_distance_matrix <- function(coords) {
  dist_matrix <- as.matrix(dist(coords))
  return(dist_matrix)
}

# Helper function to create a spatial correlation matrix
# Using the corrected inverse exponential function: Corr(d) = exp(-(d/phi_decay)^2) [cite: 101, 102]
create_correlation_matrix <- function(dist_matrix, phi_decay) {
  # Ensure phi_decay is not zero to avoid division by zero
  if (phi_decay == 0) {
    # If phi_decay is 0, implies no spatial correlation beyond identity (or error)
    # However, for positive distances, this would lead to exp(-Inf) = 0.
    # Let's return an identity matrix if phi_decay is zero, meaning independence.
    warning("phi_decay is 0, returning identity matrix for correlation.")
    return(diag(nrow(dist_matrix)))
  }
  # Corrected inverse exponential correlation function [cite: 101, 102]
  # The paper's Eq. 18 is Corr(eps_i, eps_j) = exp[-(d_ij / phi)^2]
  # where phi is the distance-decay parameter.
  corr_matrix <- exp(-(dist_matrix / phi_decay)^2)
  return(corr_matrix)
}

# Main simulation and estimation function
run_spatial_simulation_and_estimation <- function(
    N_rows = 10, N_cols = 10, # Grid dimensions
    beta_true = 1.0,          # True beta coefficient [cite: 118]
    sigma_sq_true_errors = 1.0, # True variance for generating errors (marginal variance of epsilon)
    # This corresponds to sigma^2 in the BVN assumption if Omega has 1s on diagonal
    phi_decay_errors = 1.0,   # Distance-decay parameter for errors [cite: 108]
    phi_decay_x = 1.0,        # Distance-decay parameter for X variable [cite: 115]
    sd_x_multiplier = 2.0,    # X's std dev is sd_x_multiplier * sd(errors) [cite: 113]
    num_estimation_iter = 50, # Number of iterations for BML estimation
    psi_initial_guess = 0.1   # Initial guess for psi_BML
) {
  
  n_total <- N_rows * N_cols # Total number of observations
  
  # 1. Generate coordinates
  coords <- expand.grid(x = 1:N_cols, y = 1:N_rows)
  
  # 2. Generate spatially correlated X variable
  dist_matrix_x <- calculate_distance_matrix(coords)
  R_x <- create_correlation_matrix(dist_matrix_x, phi_decay_x)
  
  # Ensure R_x is positive definite for Cholesky decomposition
  # Small jitter can help if numerically not PD
  diag(R_x) <- diag(R_x) + 1e-9 
  
  L_x <- chol(R_x)
  # Generate X such that its sd is sd_x_multiplier * sqrt(sigma_sq_true_errors)
  # Assuming the base white noise for X has variance (sd_x_multiplier * sqrt(sigma_sq_true_errors))^2
  # and correlation structure R_x
  # X_uncorr <- rnorm(n_total, mean = 0, sd = sd_x_multiplier * sqrt(sigma_sq_true_errors))
  # X_correlated <- as.vector(t(L_x) %*% X_uncorr)
  
  # Alternative: Use mvrnorm directly if R_x is a covariance matrix.
  # To make it a covariance matrix for X, scale R_x by desired variance for X
  var_x_true <- (sd_x_multiplier * sqrt(sigma_sq_true_errors))^2
  Sigma_X <- var_x_true * R_x
  # Add small value to diagonal for numerical stability if needed
  diag(Sigma_X) <- diag(Sigma_X) + 1e-9 
  X_correlated <- MASS::mvrnorm(n = 1, mu = rep(0, n_total), Sigma = Sigma_X)
  X_spatial <- as.vector(X_correlated)
  
  
  # 3. Generate spatially correlated errors (epsilon)
  dist_matrix_errors <- calculate_distance_matrix(coords)
  R_errors <- create_correlation_matrix(dist_matrix_errors, phi_decay_errors)
  
  # Ensure R_errors is positive definite
  diag(R_errors) <- diag(R_errors) + 1e-9
  
  # The simulation in paper (Eq.19) implies epsilon = Lz where z ~ N(0,1) [cite: 111, 112]
  # This means the marginal variance of epsilon_i is 1 if R_errors is a correlation matrix.
  # If sigma_sq_true_errors is meant to be the common variance factor in BVN(0, sigma^2 Omega)
  # then we scale by sqrt(sigma_sq_true_errors)
  Sigma_errors <- sigma_sq_true_errors * R_errors
  diag(Sigma_errors) <- diag(Sigma_errors) + 1e-9 # for numerical stability
  
  errors <- MASS::mvrnorm(n = 1, mu = rep(0, n_total), Sigma = Sigma_errors)
  errors_spatial <- as.vector(errors)
  
  # 4. Generate dependent variable Y
  # Y_i = beta * X_i + epsilon_i [cite: 31]
  Y_spatial <- beta_true * X_spatial + errors_spatial
  
  # 5. Center variables X and Y (as assumed in the model derivation [cite: 32])
  X_centered <- X_spatial - mean(X_spatial)
  Y_centered <- Y_spatial - mean(Y_spatial)
  
  # 6. Bivariate coding: Select pairs of neighbors
  # Strategy: horizontal pairs ( (r,c), (r,c+1) ) for odd rows and odd starting columns
  # This aims to make pairs independent of each other as per Definition 1 [cite: 39, 40]
  paired_data <- list()
  pair_idx <- 0
  for (r in 1:N_rows) {
    if (r %% 2 == 1) { # Odd rows
      for (c in 1:N_cols) {
        if (c %% 2 == 1 && (c + 1) <= N_cols) { # Odd columns, and ensure pair exists
          pair_idx <- pair_idx + 1
          
          idx1 <- (r - 1) * N_cols + c       # Index of first element in pair
          idx2 <- (r - 1) * N_cols + (c + 1) # Index of second element in pair
          
          paired_data[[pair_idx]] <- list(
            y_i = Y_centered[idx1], x_i = X_centered[idx1],
            y_l = Y_centered[idx2], x_l = X_centered[idx2]
          )
        }
      }
    }
  }
  
  q <- length(paired_data) # Number of pairs
  if (q == 0) {
    stop("No pairs were formed. Check grid dimensions and pairing strategy.")
  }
  
  # 7. Calculate sufficient statistics (alpha_1 to alpha_6) [cite: 55]
  alpha_1 <- 0; alpha_2 <- 0; alpha_3 <- 0; alpha_4 <- 0; alpha_5 <- 0; alpha_6 <- 0
  
  for (k in 1:q) {
    pair <- paired_data[[k]]
    y_i <- pair$y_i; x_i <- pair$x_i
    y_l <- pair$y_l; x_l <- pair$x_l
    
    alpha_1 <- alpha_1 + (x_i^2 + x_l^2)
    alpha_2 <- alpha_2 + (y_i^2 + y_l^2)
    alpha_3 <- alpha_3 + (x_i * y_i + x_l * y_l) # Sum of x_j*y_j for all 2q observations
    alpha_4 <- alpha_4 + (x_i * y_l + x_l * y_i) # Sum of (x_i*y_l + x_l*y_i) for each pair
    alpha_5 <- alpha_5 + (x_i * x_l)             # Sum of x_i*x_l for each pair
    alpha_6 <- alpha_6 + (y_i * y_l)             # Sum of y_i*y_l for each pair
  }
  
  # 8. Iterative BML estimation for beta, sigma_sq, psi [cite: 59]
  psi_hat <- psi_initial_guess
  beta_hat <- NA
  sigma_sq_hat <- NA
  
  cat("Starting BML estimation...\n")
  for (iter in 1:num_estimation_iter) {
    # Estimate beta_hat (Eq. 7)
    denominator_beta <- alpha_1 - 2 * psi_hat * alpha_5
    if (abs(denominator_beta) < 1e-9) {
      warning(paste("Denominator for beta_hat is near zero at iteration", iter))
      beta_hat <- ifelse(is.na(beta_hat), 0, beta_hat) # Keep last or set to 0
    } else {
      beta_hat <- (alpha_3 - psi_hat * alpha_4) / denominator_beta
    }
    
    # Estimate sigma_sq_hat (Eq. 8)
    numerator_sigma_sq <- alpha_2 + beta_hat^2 * alpha_1 - 2 * beta_hat * alpha_3 - 
      2 * psi_hat * alpha_6 - 2 * psi_hat * beta_hat^2 * alpha_5 + 
      2 * psi_hat * beta_hat * alpha_4
    
    denominator_sigma_sq <- 2 * q * (1 - psi_hat^2)
    if (abs(denominator_sigma_sq) < 1e-9 || numerator_sigma_sq < 0) {
      warning(paste("Problem with sigma_sq_hat calculation at iteration", iter,
                    " Num:", numerator_sigma_sq, "Den:", denominator_sigma_sq))
      # Keep previous sigma_sq_hat or set to a small positive if first problematic iter
      sigma_sq_hat <- ifelse(is.na(sigma_sq_hat) || sigma_sq_hat <=0, 1e-6, sigma_sq_hat) 
    } else {
      sigma_sq_hat <- numerator_sigma_sq / denominator_sigma_sq
      if (sigma_sq_hat <= 0) { # Ensure positivity
        warning(paste("sigma_sq_hat was not positive, clamping. Iter:", iter))
        sigma_sq_hat <- 1e-6 
      }
    }
    
    # Estimate psi_hat (Eq. 9)
    denominator_psi <- q * sigma_sq_hat
    if (abs(denominator_psi) < 1e-9) {
      warning(paste("Denominator for psi_hat is near zero at iteration", iter))
      # Keep previous psi_hat or default if issue from start
      psi_hat_new <- ifelse(iter==1 && psi_hat == psi_initial_guess, 0, psi_hat) 
    } else {
      psi_hat_new <- (alpha_6 - beta_hat * alpha_4 + beta_hat^2 * alpha_5) / denominator_psi
    }
    
    # Constrain psi_hat to (-1, 1)
    psi_hat_new <- max(-0.999, min(0.999, psi_hat_new))
    
    # Check for convergence (optional, for fixed iterations now)
    # if (abs(psi_hat_new - psi_hat) < 1e-6) break
    
    psi_hat <- psi_hat_new
    
    if(iter %% 10 == 0) {
      cat(paste("Iter:", iter, "beta_hat:", round(beta_hat,4), 
                "sigma_sq_hat:", round(sigma_sq_hat,4), "psi_hat:", round(psi_hat,4), "\n"))
    }
  }
  
  return(list(
    beta_hat = beta_hat,
    sigma_sq_hat = sigma_sq_hat,
    psi_hat = psi_hat,
    q_pairs = q,
    true_params = list(beta_true = beta_true, 
                       sigma_sq_true_errors = sigma_sq_true_errors, 
                       phi_decay_errors = phi_decay_errors,
                       # Calculate true psi from phi_decay_errors for adjacent cells (d=1)
                       psi_true_approx = exp(-(1/phi_decay_errors)^2) 
    )
  ))
}

# --- Example Usage ---
# Set seed for reproducibility
set.seed(123)

# Run the simulation
# The paper uses phi_decay_errors = 1 or 0.8.
# If phi_decay_errors = 1, psi for d=1 is exp(-1) = 0.367
# If phi_decay_errors = 0.8, psi for d=1 is exp(-(1/0.8)^2) = exp(-1.25^2) = exp(-1.5625) = 0.2096
# The paper reports psi = 0.367 and 0.286.
# Let's use their reported phi_decay_errors values which directly gives psi after calculation.
# For example, if target psi is 0.367 for adjacent cells (d=1), then:
# 0.367 = exp(-(1/phi_decay_errors)^2) => log(0.367) = -(1/phi_decay_errors)^2
# -log(0.367) = (1/phi_decay_errors)^2 => phi_decay_errors = 1/sqrt(-log(0.367))
# phi_decay_for_psi_0.367 <- 1 / sqrt(-log(0.367)) # approx 1.000
# phi_decay_for_psi_0.286 <- 1 / sqrt(-log(0.286)) # approx 0.882

simulation_results <- run_spatial_simulation_and_estimation(
  N_rows = 10, N_cols = 10,      # Small grid for quick test
  beta_true = 1.0,
  sigma_sq_true_errors = 1.0,   # True marginal variance of errors
  phi_decay_errors = 1.0,       # Corresponds to psi approx 0.367 for d=1
  phi_decay_x = 1.0,
  sd_x_multiplier = 2.0,
  num_estimation_iter = 20,
  psi_initial_guess = 0.2
)

cat("\n--- Simulation Results ---\n")
print(simulation_results)

# To run for parameters similar to Experiment 1 in paper [cite: 122] (n=100, psi=0.367):
# results_exp1_like <- run_spatial_simulation_and_estimation(
#   N_rows = 10, N_cols = 10,
#   beta_true = 1.0,
#   sigma_sq_true_errors = 1.0, # Assuming true error variance for BVN is 1
#   phi_decay_errors = 1/sqrt(-log(0.367)), # To get psi ~ 0.367 for d=1
#   phi_decay_x = 1.0, # Paper states phi=1 for x in all experiments [cite: 115]
#   sd_x_multiplier = 2.0,
#   num_estimation_iter = 100, # More iterations for stability
#   psi_initial_guess = 0.3
# )
# print(results_exp1_like)
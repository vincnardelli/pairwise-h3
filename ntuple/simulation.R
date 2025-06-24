# Install and load necessary packages
# Check if packages are installed, if not, install them
if (!requireNamespace("spdep", quietly = TRUE)) install.packages("spdep")
if (!requireNamespace("nloptr", quietly = TRUE)) install.packages("nloptr") # For numerical optimization
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")
if (!requireNamespace("combinat", quietly = TRUE)) install.packages("combinat") # For combn (ensure it's available)
if (!requireNamespace("here", quietly = TRUE)) install.packages("here") # For path management
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr") # For read_rds/write_rds

# Load required libraries
library(spdep)    # For spatial weights and SAR model utilities
library(nloptr)   # For numerical optimization (e.g., L-BFGS-B algorithm)
library(ggplot2)  # For plotting results
library(reshape2) # For data manipulation (e.g., melt for ggplot)
library(combinat) # For combn
library(here)     # For convenient path management
library(readr)    # For saving/loading R objects

# --- Utility Function: Generate Spatial Weights Matrix (QUEEN Contiguity - CORRECTED) ---
# This function creates a spatial weights matrix based on a grid.
# Changed to "queen" contiguity to allow for 3-cliques (triangles).
generate_grid_W <- function(n_rows, n_cols) {
  N <- n_rows * n_cols
  coords <- expand.grid(x = 1:n_cols, y = 1:n_rows)
  lw <- cell2nb(n_cols, n_rows, type = "queen") # CHANGED FROM "rook" to "queen" 
  W <- nb2mat(lw, style = "B", zero.policy = TRUE) # Convert to binary matrix
  return(W)
}

# --- Utility Function: Simulate SAR Data ---
# Generates data according to a SAR process: y = (I - rho*W)^-1 * (X*beta + epsilon)
simulate_sar_data <- function(rho, beta, sigma2, W, X) {
  N <- nrow(W)
  I_N <- diag(N)
  epsilon <- rnorm(N, mean = 0, sd = sqrt(sigma2))
  # Inverse term: (I - rho*W)^-1
  inv_term <- solve(I_N - rho * W)
  y <- inv_term %*% (X %*% beta + epsilon)
  return(list(y = as.vector(y), X = X, W = W))
}

# --- Utility Function: Create Clusters (CORRECTED for n=3 clique formation) ---
# This function creates overlapping clusters, mimicking real-world composite likelihood applications.
# For n=3, it now specifically identifies 3-cliques (triangles) in the spatial graph.
create_clusters <- function(N, n_cluster_size, W) {
  clusters <- list()
  # Using a hash-set-like structure to store unique sorted clusters for efficiency
  unique_clusters_set <- new.env(hash = TRUE, parent = emptyenv())
  
  if (n_cluster_size == 2) {
    for (i in 1:N) {
      neighbors <- which(W[i, ] == 1)
      for (j in neighbors) {
        if (i < j) { # Avoid duplicates (i,j) and (j,i)
          sorted_pair <- paste(sort(c(i, j)), collapse = "_")
          if (is.null(unique_clusters_set[[sorted_pair]])) {
            clusters[[length(clusters) + 1]] <- c(i, j)
            unique_clusters_set[[sorted_pair]] <- TRUE
          }
        }
      }
    }
  } else if (n_cluster_size == 3) {
    # For triplets, ensure (i, j, k) form a clique based on W (i.e., each is a neighbor of the other two)
    # This is critical for consistency with the Jacobian and epsilon definition in the paper.
    for (i in 1:N) {
      neighbors_i <- which(W[i, ] == 1)
      if (length(neighbors_i) >= 2) {
        # Iterate over all distinct pairs of neighbors of i
        for (j_idx in 1:(length(neighbors_i)-1)) {
          j <- neighbors_i[j_idx]
          for (k_idx in (j_idx+1):length(neighbors_i)) {
            k <- neighbors_i[k_idx]
            # Crucially, check if j and k are also neighbors of each other
            if (W[j, k] == 1) { # This condition ensures (i,j,k) forms a clique 
              sorted_triplet <- paste(sort(c(i, j, k)), collapse = "_")
              if (is.null(unique_clusters_set[[sorted_triplet]])) {
                clusters[[length(clusters) + 1]] <- c(i, j, k)
                unique_clusters_set[[sorted_triplet]] <- TRUE
              }
            }
          }
        }
      }
    }
  }
  return(clusters)
}

# --- CML Estimator for n=2 (Pairwise) ---
cml_n2_estimator <- function(y, X, clusters, N_total) {
  q <- length(clusters) # Number of pairs
  
  # Calculate sufficient statistics
  T1 = 0; T2 = 0; T3 = 0; T4 = 0; T5 = 0
  for (k in 1:q) {
    i_k <- clusters[[k]][1]
    l_k <- clusters[[k]][2]
    T1 <- T1 + (X[i_k]^2 + X[l_k]^2) # Assuming X is a vector for a single regressor
    T2 <- T2 + (y[i_k]^2 + y[l_k]^2)
    T3 <- T3 + (X[i_k]*y[i_k] + X[l_k]*y[l_k])
    T4 <- T4 + (y[i_k]*y[l_k])
    T5 <- T5 + (X[i_k]*y[l_k] + X[l_k]*y[i_k])
  }
  
  # Log-likelihood function for optimization (concentrated on rho)
  loglik_n2 <- function(rho_param) {
    rho <- rho_param[1]
    
    # Jacobian term (from your paper)
    jac_term <- (1 - rho^2)
    if (jac_term <= 0) return(-Inf) # Avoid log of non-positive numbers
    
    # Conditional beta estimator
    beta_hat <- (T3 - rho * T5) / T1
    
    # Conditional sigma^2 estimator
    sigma2_hat <- (1 / (2 * q)) * (
      (1 + rho^2) * T2 - 4 * rho * T4 - 2 * beta_hat * T3 + 2 * rho * beta_hat * T5 + beta_hat^2 * T1
    )
    if (sigma2_hat <= 0) return(-Inf) # sigma^2 must be positive
    
    # Total log-likelihood
    val <- -(2 * q / 2) * log(2 * pi * sigma2_hat) + q * log(abs(jac_term)) - (1 / (2 * sigma2_hat)) * (2 * q * sigma2_hat)
    return(-val) # nloptr minimizes, so return negative log-likelihood
  }
  
  # Initial values for optimization
  initial_rho <- 0.1 # Sensible starting point
  
  # Optimize using L-BFGS-B
  opts <- list("algorithm" = "NLOPT_LN_SBPLX", "xtol_rel" = 1.0e-8, "maxeval" = 1000)
  result <- nloptr(x0 = initial_rho, eval_f = loglik_n2,
                   lb = -0.99, ub = 0.99, # Bounds for rho
                   opts = opts)
  
  rho_hat <- result$solution[1]
  beta_hat <- (T3 - rho_hat * T5) / T1
  sigma2_hat <- (1 / (2 * q)) * (
    (1 + rho_hat^2) * T2 - 4 * rho_hat * T4 - 2 * beta_hat * T3 + 2 * rho_hat * beta_hat * T5 + beta_hat^2 * T1
  )
  
  return(list(rho = rho_hat, beta = beta_hat, sigma2 = sigma2_hat, converged = result$status))
}

# --- CML Estimator for n=3 (Triplet-wise) ---
cml_n3_estimator <- function(y, X, clusters, N_total) {
  q <- length(clusters) # Number of triplets
  
  # Calculate sufficient statistics (S1 to S9)
  S1 = 0; S2 = 0; S3 = 0; S4 = 0; S5 = 0; S6 = 0; S7 = 0; S8 = 0; S9 = 0
  for (k in 1:q) {
    i_k <- clusters[[k]][1]
    l_k <- clusters[[k]][2]
    m_k <- clusters[[k]][3]
    
    S1 <- S1 + (X[i_k]^2 + X[l_k]^2 + X[m_k]^2)
    S2 <- S2 + (y[i_k]^2 + y[l_k]^2 + y[m_k]^2)
    S3 <- S3 + (X[i_k]*y[i_k] + X[l_k]*y[l_k] + X[m_k]*y[m_k])
    S4 <- S4 + (y[i_k]*y[l_k])
    S5 <- S5 + (y[i_k]*y[m_k])
    S6 <- S6 + (y[l_k]*y[m_k])
    S7 <- S7 + (X[i_k]*y[l_k] + X[l_k]*y[i_k])
    S8 <- S8 + (X[i_k]*y[m_k] + X[m_k]*y[i_k])
    S9 <- S9 + (X[l_k]*y[m_k] + X[m_k]*y[l_k])
  }
  
  # Log-likelihood function for optimization (concentrated on rho)
  loglik_n3 <- function(rho_param) {
    rho <- rho_param[1]
    
    # Jacobian term (from your paper)
    # det(J_3*) = 1 - 3*rho^2 - 2*rho^3
    jac_term <- (1 - 3*rho^2 - 2*rho^3)
    if (jac_term <= 0) return(-Inf) # Avoid log of non-positive numbers
    
    # Conditional beta estimator
    beta_hat <- (S3 - rho * (S7 + S8 + S9)) / S1
    
    # Conditional sigma^2 estimator (CORRECTED)
    sigma2_hat <- (1 / (3 * q)) * (
      (1 + 2*rho^2)*S2 + (-4*rho + 2*rho^2)*(S4+S5+S6) - # CORRECTED TERM 
        2*beta_hat*S3 + 2*rho*beta_hat*(S7+S8+S9) + beta_hat^2*S1
    )
    if (sigma2_hat <= 0) return(-Inf) # sigma^2 must be positive
    
    # Total log-likelihood
    val <- -(3 * q / 2) * log(2 * pi * sigma2_hat) + q * log(abs(jac_term)) - (1 / (2 * sigma2_hat)) * (3 * q * sigma2_hat)
    return(-val) # nloptr minimizes
  }
  
  # Initial values for optimization
  initial_rho <- 0.1
  
  # Optimize using L-BFGS-B
  # Bounds for rho: Jacobian must be positive.
  # 1 - 3*rho^2 - 2*rho^3 > 0. Roots are approx -1.618, -0.618, 0.5.
  # So, valid range within the typical SAR (-1,1) is approximately (-0.618, 0.5)
  # Let's use (-0.6, 0.49) as a conservative bound from your paper's Jacobian for n=3.
  opts <- list("algorithm" = "NLOPT_LN_SBPLX", "xtol_rel" = 1.0e-8, "maxeval" = 1000)
  result <- nloptr(x0 = initial_rho, eval_f = loglik_n3,
                   lb = -0.6, ub = 0.49, # Conservative bounds for rho for n=3 Jacobian
                   opts = opts)
  
  rho_hat <- result$solution[1]
  beta_hat <- (S3 - rho_hat * (S7 + S8 + S9)) / S1
  # Final sigma2_hat calculation outside loglik_n3, using the optimized rho_hat
  sigma2_hat <- (1 / (3 * q)) * (
    (1 + 2*rho_hat^2)*S2 + (-4*rho_hat + 2*rho_hat^2)*(S4+S5+S6) - # CORRECTED TERM
      2*beta_hat*S3 + 2*rho_hat*beta_hat*(S7+S8+S9) + beta_hat^2*S1
  )
  
  return(list(rho = rho_hat, beta = beta_hat, sigma2 = sigma2_hat, converged = result$status))
}

# --- Simulation Parameters ---
N_rows <- 20 # Number of rows in the grid
N_cols <- 20 # Number of columns in the grid
N_total <- N_rows * N_cols # Total number of spatial units (e.g., 400)

True_rho <- 0.4    # True spatial autoregressive parameter
True_beta <- 2.5   # True regression coefficient
True_sigma2 <- 1.0 # True error variance

Num_replications <- 100 # Number of Monte Carlo simulations (increase for more robust results, e.g., 500-1000)

# Prepare spatial weights matrix once (now using QUEEN contiguity)
W_matrix <- generate_grid_W(N_rows, N_cols)

# Prepare covariates (single regressor X)
set.seed(123) # For reproducibility
X_matrix <- matrix(rnorm(N_total, mean = 10, sd = 2), ncol = 1) # Single regressor

# Create directories for saving results if they don't exist
# This is assuming your working directory is the root of your project
# and 'ntuple/data' and 'ntuple/images' are subfolders.
if (!dir.exists(here("ntuple", "data"))) {
  dir.create(here("ntuple", "data"), recursive = TRUE)
}
if (!dir.exists(here("ntuple", "images"))) {
  dir.create(here("ntuple", "images"), recursive = TRUE)
}


# Store results
results_n2 <- data.frame(rho_hat = numeric(Num_replications),
                         beta_hat = numeric(Num_replications),
                         sigma2_hat = numeric(Num_replications),
                         time_taken = numeric(Num_replications))

results_n3 <- data.frame(rho_hat = numeric(Num_replications),
                         beta_hat = numeric(Num_replications),
                         sigma2_hat = numeric(Num_replications),
                         time_taken = numeric(Num_replications))

# --- Main Simulation Loop ---
cat("Starting simulation...\n")
for (rep in 1:Num_replications) {
  cat(sprintf("Replication %d/%d\n", rep, Num_replications))
  
  # Generate SAR data for this replication
  sim_data <- simulate_sar_data(True_rho, True_beta, True_sigma2, W_matrix, X_matrix)
  y_sim <- sim_data$y
  X_sim <- sim_data$X
  
  # --- CML n=2 (Pairwise) Estimation ---
  start_time_n2 <- Sys.time()
  clusters_n2 <- create_clusters(N_total, 2, W_matrix)
  est_n2 <- cml_n2_estimator(y_sim, X_sim, clusters_n2, N_total)
  end_time_n2 <- Sys.time()
  results_n2[rep, "time_taken"] <- as.numeric(difftime(end_time_n2, start_time_n2, units = "secs"))
  results_n2[rep, "rho_hat"] <- est_n2$rho
  results_n2[rep, "beta_hat"] <- est_n2$beta
  results_n2[rep, "sigma2_hat"] <- est_n2$sigma2
  
  # --- CML n=3 (Triplet-wise) Estimation ---
  start_time_n3 <- Sys.time()
  clusters_n3 <- create_clusters(N_total, 3, W_matrix) # Will now select cliques only
  # Only proceed if clusters_n3 is not empty (it will be if no 3-cliques are found)
  if (length(clusters_n3) > 0) {
    est_n3 <- cml_n3_estimator(y_sim, X_sim, clusters_n3, N_total)
    end_time_n3 <- Sys.time()
    results_n3[rep, "time_taken"] <- as.numeric(difftime(end_time_n3, start_time_n3, units = "secs"))
    results_n3[rep, "rho_hat"] <- est_n3$rho
    results_n3[rep, "beta_hat"] <- est_n3$beta
    results_n3[rep, "sigma2_hat"] <- est_n3$sigma2
  } else {
    warning(paste("No triplets (cliques) formed for N =", N_total, ". Skipping n=3 estimation for replication", rep))
    results_n3[rep, c("rho_hat", "beta_hat", "sigma2_hat", "time_taken")] <- NA
  }
}
cat("Simulation finished!\n")

# --- Save results for later analysis ---
# write_rds(results_n2, here("ntuple", "data", "results_n2.rds"))
# write_rds(results_n3, here("ntuple", "data", "results_n3.rds"))
cat(paste0("Simulation results saved to ", here("ntuple", "data"), "\n"))


# --- Results Analysis ---
calculate_summary <- function(results_df, true_val, param_name) {
  estimates <- results_df[[paste0(param_name, "_hat")]]
  estimates <- estimates[!is.na(estimates)] # Remove NAs if any (e.g., from n=3 skipping)
  
  if (length(estimates) == 0) {
    return(data.frame(
      Parameter = param_name,
      TrueValue = true_val,
      MeanEstimate = NA,
      Bias = NA,
      Variance = NA,
      MSE = NA
    ))
  }
  
  bias <- mean(estimates) - true_val
  variance <- var(estimates)
  mse <- mean((estimates - true_val)^2)
  
  return(data.frame(
    Parameter = param_name,
    TrueValue = true_val,
    MeanEstimate = mean(estimates),
    Bias = bias,
    Variance = variance,
    MSE = mse
  ))
}

# --- Load results ANALYSIS -----
# If you run the script continuously, results_n2 and results_n3 are already in memory.
# If you restart R, uncomment these lines and ensure 'here' package points to your project root.
results_n2 <- read_rds(here("ntuple", "data", "results_n2.rds"))
results_n3 <- read_rds(here("ntuple", "data", "results_n3.rds"))

# Summarize results for n=2
summary_n2_rho <- calculate_summary(results_n2, True_rho, "rho")
summary_n2_beta <- calculate_summary(results_n2, True_beta, "beta")
summary_n2_sigma2 <- calculate_summary(results_n2, True_sigma2, "sigma2")

# Summarize results for n=3 (handling NAs from failed estimations, if any)
results_n3_clean <- results_n3[complete.cases(results_n3), ]
summary_n3_rho <- calculate_summary(results_n3_clean, True_rho, "rho")
summary_n3_beta <- calculate_summary(results_n3_clean, True_beta, "beta")
summary_n3_sigma2 <- calculate_summary(results_n3_clean, True_sigma2, "sigma2")

# Combine summaries
summary_table <- rbind(
  cbind(Model = "Pairwise (n=2)", summary_n2_rho),
  cbind(Model = "Pairwise (n=2)", summary_n2_beta),
  cbind(Model = "Pairwise (n=2)", summary_n2_sigma2),
  cbind(Model = "Triplet-wise (n=3)", summary_n3_rho),
  cbind(Model = "Triplet-wise (n=3)", summary_n3_beta),
  cbind(Model = "Triplet-wise (n=3)", summary_n3_sigma2)
)

print("--- Parameter Estimation Summary ---")
print(summary_table)

# Compare computational times
mean_time_n2 <- mean(results_n2$time_taken)
mean_time_n3 <- mean(results_n3$time_taken, na.rm = TRUE)

cat("\n--- Computational Time Comparison ---\n")
cat(sprintf("Mean time for Pairwise (n=2) estimation: %.4f seconds\n", mean_time_n2))
cat(sprintf("Mean time for Triplet-wise (n=3) estimation: %.4f seconds\n", mean_time_n3))

# --- Visualization ---

# Reshape data for plotting
plot_data_rho <- data.frame(
  Replication = 1:Num_replications,
  Pairwise = results_n2$rho_hat,
  Tripletwise = results_n3$rho_hat
)
plot_data_rho_melted <- melt(plot_data_rho, id.vars = "Replication", variable.name = "Estimator", value.name = "Rho_Estimate")

plot_data_beta <- data.frame(
  Replication = 1:Num_replications,
  Pairwise = results_n2$beta_hat,
  Tripletwise = results_n3$beta_hat
)
plot_data_beta_melted <- melt(plot_data_beta, id.vars = "Replication", variable.name = "Estimator", value.name = "Beta_Estimate")

plot_data_sigma2 <- data.frame(
  Replication = 1:Num_replications,
  Pairwise = results_n2$sigma2_hat,
  Tripletwise = results_n3$sigma2_hat
)
plot_data_sigma2_melted <- melt(plot_data_sigma2, id.vars = "Replication", variable.name = "Estimator", value.name = "Sigma2_Estimate")


# Plotting rho estimates
ggplot(plot_data_rho_melted, aes(x = Estimator, y = Rho_Estimate, fill = Estimator)) +
  geom_boxplot() +
  geom_hline(yintercept = True_rho, linetype = "dashed", color = "red", linewidth = 1) + # Changed size to linewidth for modern ggplot
  labs(title = bquote("Distribution of " ~ rho ~ " Estimates (" ~ .(Num_replications) ~ " Replications)"),
       y = expression(rho ~ "Estimate"),
       x = "CML Estimator Type") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")
ggsave(here("ntuple","images","plot_data_rho_melted.pdf"), width = 8, height = 5) # Adjusted size for better fit in paper


# Plotting beta estimates
ggplot(plot_data_beta_melted, aes(x = Estimator, y = Beta_Estimate, fill = Estimator)) +
  geom_boxplot() +
  geom_hline(yintercept = True_beta, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = bquote("Distribution of " ~ beta ~ " Estimates (" ~ .(Num_replications) ~ " Replications)"),
       y = expression(beta ~ "Estimate"),
       x = "CML Estimator Type") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")
ggsave(here("ntuple","images","plot_data_beta_melted.pdf"), width = 8, height = 5)


# Plotting sigma2 estimates
ggplot(plot_data_sigma2_melted, aes(x = Estimator, y = Sigma2_Estimate, fill = Estimator)) +
  geom_boxplot() +
  geom_hline(yintercept = True_sigma2, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = bquote("Distribution of " ~ sigma^2 ~ " Estimates (" ~ .(Num_replications) ~ " Replications)"),
       y = expression(sigma^2 ~ "Estimate"),
       x = "CML Estimator Type") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")
ggsave(here("ntuple","images","plot_data_sigma2_melted.pdf"), width = 8, height = 5)


# Plotting computational times
time_data <- data.frame(
  Model = c("Pairwise (n=2)", "Triplet-wise (n=3)"),
  MeanTime = c(mean_time_n2, mean_time_n3)
)

ggplot(time_data, aes(x = Model, y = MeanTime, fill = Model)) +
  geom_bar(stat = "identity") +
  labs(title = "Mean Computational Time per Replication",
       y = "Time (seconds)",
       x = "CML Estimator Type") +
  theme_minimal() +
  scale_fill_brewer(palette = "Pastel1")
ggsave(here("ntuple","images","time_data.pdf"), width = 8, height = 5)

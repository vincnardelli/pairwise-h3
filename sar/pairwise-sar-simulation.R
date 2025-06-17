#----------------------------------------------------------------------#
# SCRIPT DI SIMULAZIONE PER PAIRWISE LIKELIHOOD SAR (v2 - CORRETTO)
# Implementa il sub-campionamento con "coding" di Besag/Arbia
#----------------------------------------------------------------------#
library(spdep)
library(spatialreg)
library(Matrix)
library(dplyr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(here)
library(readr)

# --- 2. Parametri della Simulazione ----
set.seed(42) # Per la riproducibilità

# Parametri della griglia
GRID_SIZE <- 25       # Griglia 25x25 -> N=625
N <- GRID_SIZE^2

# Parametri veri del modello
TRUE_BETA0 <- 1.0
TRUE_BETA1 <- 2.0
TRUE_RHO   <- 0.7
TRUE_SIGMA2<- 1.5

# Parametri della simulazione
N_SIM <- 500          # Numero di ripetizioni Monte Carlo

# --- 3. Funzione per Generare Dati SAR su Griglia (invariata) ----
generate_sar_data_grid <- function(grid_size, beta_true, rho_true, sigma2_true) {
  n_points <- grid_size^2
  nb <- spdep::cell2nb(grid_size, grid_size, type = "queen")
  W_listw <- spdep::nb2listw(nb, style = "W")
  W <- as(W_listw, "CsparseMatrix")
  
  x <- rnorm(n_points, mean = 10, sd = 3)
  X_matrix <- cbind(1, x)
  epsilon <- rnorm(n_points, mean = 0, sd = sqrt(sigma2_true))
  
  I <- .sparseDiagonal(n_points)
  A <- I - rho_true * W
  y <- solve(A, X_matrix %*% beta_true + epsilon)
  
  return(list(y = as.vector(y), x = x, nb_list = nb, listw = W_listw))
}

# --- 4. NUOVA Funzione per creare coppie non sovrapposte (BIVARIATE CODING) ----
create_coding_pairs <- function(grid_size, nb_list, start_row = 1, start_col = 1) {
  
  coords <- expand.grid(row = 1:grid_size, col = 1:grid_size)
  coords$id <- 1:(grid_size^2)
  coords$is_used <- FALSE
  
  selected_pairs <- list()
  
  # Itera sulla griglia seguendo un pattern a scacchiera
  for (r in seq(start_row, grid_size, by = 2)) {
    for (c in seq(start_col, grid_size, by = 2)) {
      
      origin_index <- coords$id[coords$row == r & coords$col == c]
      
      # Se l'origine è già stata usata (come vicino di un'altra coppia), salta
      if (coords$is_used[origin_index]) next
      
      # Trova i vicini non ancora utilizzati
      potential_neighbors <- nb_list[[origin_index]]
      available_neighbors <- potential_neighbors[!coords$is_used[potential_neighbors]]
      
      if (length(available_neighbors) > 0) {
        # Scegli un vicino a caso tra quelli disponibili
        neighbor_index <- sample(available_neighbors, 1)
        
        # Aggiungi la coppia alla lista
        selected_pairs[[length(selected_pairs) + 1]] <- data.frame(i = origin_index, l = neighbor_index)
        
        # Marca entrambe le osservazioni come usate
        coords$is_used[origin_index] <- TRUE
        coords$is_used[neighbor_index] <- TRUE
      }
    }
  }
  
  if (length(selected_pairs) > 0) {
    return(do.call(rbind, selected_pairs))
  } else {
    return(data.frame(i = integer(), l = integer()))
  }
}

# --- 5. Funzione di Stima (MODIFICATA per accettare le coppie) ----
estimate_pairwise_sar <- function(y, x, pairs_df) {
  
  q <- nrow(pairs_df)
  if (q == 0) return(c(beta_hat = NA, rho_hat = NA, sigma2_hat = NA))
  
  y_c <- y - mean(y); x_c <- x - mean(x)
  
  xi <- x_c[pairs_df$i]; xl <- x_c[pairs_df$l]
  yi <- y_c[pairs_df$i]; yl <- y_c[pairs_df$l]
  
  T1 <- sum(xi^2 + xl^2); T2 <- sum(yi^2 + yl^2); T3 <- sum(xi*yi + xl*yl)
  T4 <- sum(yi*yl); T5 <- sum(xi*yl + xl*yi)
  
  rho_grid <- seq(-0.99, 0.99, by = 0.01)
  log_lik_values <- numeric(length(rho_grid))
  
  for (i in 1:length(rho_grid)) {
    rho <- rho_grid[i]
    beta <- (T3 - rho * T5) / T1
    sigma2 <- (1/(2*q)) * ((1+rho^2)*T2 + beta^2*T1 - 2*beta*T3 - 4*rho*T4 + 2*rho*beta*T5)
    
    if (is.na(sigma2) || sigma2 <= 0 || (1 - rho^2) <= 0) {
      log_lik_values[i] <- -Inf
    } else {
      log_lik_values[i] <- -q * log(sigma2) + q * log(1 - rho^2)
    }
  }
  
  best_rho <- rho_grid[which.max(log_lik_values)]
  final_beta <- (T3 - best_rho * T5) / T1
  final_sigma2 <- (1/(2*q)) * ((1+best_rho^2)*T2 + final_beta^2*T1 - 2*final_beta*T3 - 4*best_rho*T4 + 2*best_rho*final_beta*T5)
  
  return(c(beta_hat = final_beta, rho_hat = best_rho, sigma2_hat = final_sigma2))
}


# --- 6. Loop di Simulazione  ----
results_df <- data.frame(
  beta_pw = numeric(N_SIM), rho_pw = numeric(N_SIM), sigma2_pw = numeric(N_SIM), time_pw = numeric(N_SIM),
  beta_ml = numeric(N_SIM), rho_ml = numeric(N_SIM), sigma2_ml = numeric(N_SIM), time_ml = numeric(N_SIM)
)

cat(paste("Avvio simulazione (con CODING CORRETTO) per N =", N, "e", N_SIM, "ripetizioni...\n"))
progress_bar <- txtProgressBar(min = 0, max = N_SIM, style = 3)

for (i in 1:N_SIM) {
  sim_data <- generate_sar_data_grid(GRID_SIZE, c(TRUE_BETA0, TRUE_BETA1), TRUE_RHO, TRUE_SIGMA2)
  
  # --- Stima con Pairwise SAR (con media su 4 schemi di coding) ---
  start_time_pw <- Sys.time()
  
  # I 4 possibili schemi di coding a scacchiera (come da Arbia, 2014)
  coding_schemes <- list(
    create_coding_pairs(GRID_SIZE, sim_data$nb_list, start_row = 1, start_col = 1),
    create_coding_pairs(GRID_SIZE, sim_data$nb_list, start_row = 1, start_col = 2),
    create_coding_pairs(GRID_SIZE, sim_data$nb_list, start_row = 2, start_col = 1),
    create_coding_pairs(GRID_SIZE, sim_data$nb_list, start_row = 2, start_col = 2)
  )
  
  # Esegui la stima per ogni schema e fai la media
  estimates_from_schemes <- lapply(coding_schemes, function(pairs) {
    estimate_pairwise_sar(sim_data$y, sim_data$x, pairs)
  })
  
  # Calcola la media delle stime valide
  valid_estimates <- do.call(rbind, estimates_from_schemes)
  avg_pw_estimates <- colMeans(valid_estimates, na.rm = TRUE)
  
  end_time_pw <- Sys.time()
  
  results_df$beta_pw[i] <- avg_pw_estimates["beta_hat"]
  results_df$rho_pw[i] <- avg_pw_estimates["rho_hat"]
  results_df$sigma2_pw[i] <- avg_pw_estimates["sigma2_hat"]
  results_df$time_pw[i] <- as.numeric(end_time_pw - start_time_pw)
  
  # --- Stima con ML standard (benchmark) ---
  start_time_ml <- Sys.time()
  df_for_ml <- data.frame(y = sim_data$y, x = sim_data$x)
  ml_model <- spatialreg::lagsarlm(y ~ x, data = df_for_ml, listw = sim_data$listw, type = "lag", zero.policy=TRUE)
  end_time_ml <- Sys.time()
  
  results_df$beta_ml[i] <- coef(ml_model)["x"]
  results_df$rho_ml[i] <- ml_model$rho
  results_df$sigma2_ml[i] <- ml_model$s2
  results_df$time_ml[i] <- as.numeric(end_time_ml - start_time_ml)
  
  setTxtProgressBar(progress_bar, i)
}
close(progress_bar)
cat("\nSimulazione completata.\n\n")

# readr::write_rds(x = results_df, file = here("sar", "data", "results_df.rds"))
results_df = readr::read_rds(file =  here("sar", "data", "results_df.rds"))


# --- 7. Analisi e Stampa dei Risultati (invariato) ----
results_clean <- na.omit(results_df)
bias_beta_pw <- mean(results_clean$beta_pw) - TRUE_BETA1
mse_beta_pw <- mean((results_clean$beta_pw - TRUE_BETA1)^2)
bias_rho_pw <- mean(results_clean$rho_pw) - TRUE_RHO
mse_rho_pw <- mean((results_clean$rho_pw - TRUE_RHO)^2)

bias_beta_ml <- mean(results_clean$beta_ml) - TRUE_BETA1
mse_beta_ml <- mean((results_clean$beta_ml - TRUE_BETA1)^2)
bias_rho_ml <- mean(results_clean$rho_ml) - TRUE_RHO
mse_rho_ml <- mean((results_clean$rho_ml - TRUE_RHO)^2)

summary_table <- data.frame(
  Parameter = c("Beta1", "Rho", "Avg Time (s)"),
  True_Value = c(TRUE_BETA1, TRUE_RHO, NA),
  Pairwise_SAR_Mean = c(mean(results_clean$beta_pw), mean(results_clean$rho_pw), mean(results_clean$time_pw)),
  Pairwise_SAR_Bias = c(bias_beta_pw, bias_rho_pw, NA),
  Pairwise_SAR_MSE = c(mse_beta_pw, mse_rho_pw, NA),
  Standard_ML_Mean = c(mean(results_clean$beta_ml), mean(results_clean$rho_ml), mean(results_clean$time_ml)),
  Standard_ML_Bias = c(bias_beta_ml, bias_rho_ml, NA),
  Standard_ML_MSE = c(mse_beta_ml, mse_rho_ml, NA)
)



# --- 7. Grafici ----
## --- 7.2 Preparazione dei Dati per ggplot2 ----
# Trasformiamo i dati da un formato "wide" a "long" per plottarli facilmente

# Dati per i parametri Beta e Rho
results_long_params <- results_clean %>%
  select(beta_pw, rho_pw, beta_ml, rho_ml) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("parameter", "method"),
    names_pattern = "(beta|rho)_(pw|ml)",
    values_to = "estimate"
  ) %>%
  mutate(
    method = factor(method, levels = c("pw", "ml"), labels = c("Pairwise SAR", "Standard ML")),
    parameter = factor(parameter, levels = c("beta", "rho"), labels = c("Beta (Slope)", "Rho (Spatial Corr.)"))
  )

# Aggiungiamo i valori veri per i plot
true_values <- data.frame(
  parameter = factor(c("Beta (Slope)", "Rho (Spatial Corr.)"), levels = c("Beta (Slope)", "Rho (Spatial Corr.)")),
  value = c(TRUE_BETA1, TRUE_RHO)
)

# Dati per il tempo di calcolo
results_long_time <- results_clean %>%
  select(time_pw, time_ml) %>%
  pivot_longer(
    cols = everything(),
    names_to = "method",
    names_pattern = "time_(pw|ml)",
    values_to = "time_seconds"
  ) %>%
  mutate(
    method = factor(method, levels = c("pw", "ml"), labels = c("Pairwise SAR", "Standard ML"))
  )

## --- 7.3. Creazione dei Grafici ----

# Grafico 1: Boxplot delle distribuzioni degli stimatori
p_estimates <- ggplot(results_long_params, aes(x = method, y = estimate, fill = method)) +
  geom_boxplot(show.legend = FALSE, outlier.shape = 21, outlier.size = 1.5, outlier.alpha = 0.3) +
  geom_hline(data = true_values, aes(yintercept = value), linetype = "dashed", color = "red", size = 1) +
  facet_wrap(~parameter, scales = "free_y") +
  labs(
    title = "Distribution of Parameter Estimates",
    subtitle = "Comparison of Pairwise SAR vs. Standard Maximum Likelihood",
    x = "Estimation Method",
    y = "Estimated Value"
  ) +
  scale_fill_brewer(palette = "Pastel1") +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )

# Grafico 2: Boxplot dei tempi di calcolo (con scala logaritmica)
p_time <- ggplot(results_long_time, aes(x = method, y = time_seconds, fill = method)) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", show.legend = FALSE, outlier.alpha = 0.5) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(
    title = "Computational Time",
    subtitle = "Note: Y-axis is on a logarithmic scale",
    x = "Estimation Method",
    y = "Time (seconds)"
  ) +
  scale_fill_brewer(palette = "Pastel1") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 11, face = "italic")
  )

# Grafico 3: Correlazione tra gli stimatori Pairwise (per confermare la teoria)
p_correlation <- ggplot(results_clean, aes(x = beta_pw, y = rho_pw)) +
  geom_point(alpha = 0.2, color = "darkblue") +
  geom_density_2d(color = "black") +
  geom_smooth(method = "lm", color = "firebrick", se = FALSE) +
  labs(
    title = "Estimator Correlation in Pairwise SAR",
    subtitle = "Confirms the non-diagonal Fisher Information Matrix",
    x = "Beta (Slope) Estimate",
    y = "Rho Estimate"
  ) +
  annotate("text", x = min(results_clean$beta_pw), y = max(results_clean$rho_pw), 
           label = paste("Correlation:", round(cor(results_clean$beta_pw, results_clean$rho_pw), 2)),
           hjust = 0, vjust = 1, size = 5, fontface="bold") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )

# Visualizza i grafici
p_estimates
p_time
p_correlation

# Combinare i grafici in un unico layout e salvarli
# (il layout può essere modificato a piacere)
final_plot_layout <- (p_estimates | p_time)

# Salva il layout combinato
ggsave(here("sar","images","p_estimates.pdf"), plot = p_estimates, width = 16, height = 8)
ggsave(here("sar","images","p_time.pdf"), plot = p_time, width = 16, height = 8)
ggsave(here("sar","images","simulation_summary_plot.pdf"), plot = final_plot_layout, width = 16, height = 8)
ggsave(here("sar","images","correlation_plot.pdf"), plot = p_correlation, width = 8, height = 7)

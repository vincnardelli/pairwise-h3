# Carica le librerie necessarie
# install.packages("MASS") # Se non già installato
library(MASS) # Per mvrnorm

#-------------------------------------------------------------------------------
# Funzioni Helper (dallo script precedente)
#-------------------------------------------------------------------------------

# Funzione helper per calcolare la matrice delle distanze Euclidee
calculate_distance_matrix <- function(coords) {
  dist_matrix <- as.matrix(dist(coords))
  return(dist_matrix)
}

# Funzione helper per creare una matrice di correlazione spaziale
# Utilizzando la funzione esponenziale inversa corretta: Corr(d) = exp(-(d/phi_decay)^2)
create_correlation_matrix <- function(dist_matrix, phi_decay) {
  if (phi_decay == 0) {
    warning("phi_decay è 0, restituisce una matrice identità per la correlazione.")
    return(diag(nrow(dist_matrix)))
  }
  corr_matrix <- exp(-(dist_matrix / phi_decay)^2)
  return(corr_matrix)
}

#-------------------------------------------------------------------------------
# Funzione di Simulazione e Stima Singola (modificata per output più pulito)
#-------------------------------------------------------------------------------
run_spatial_simulation_and_estimation <- function(
    N_rows = 10, N_cols = 10,    # Dimensioni della griglia
    beta_true = 1.0,            # Valore vero di beta
    sigma_sq_true_errors = 1.0, # Varianza vera per la generazione degli errori (varianza marginale di epsilon)
    phi_decay_errors = 1.0,     # Parametro di decadimento della distanza per gli errori
    phi_decay_x = 1.0,          # Parametro di decadimento della distanza per la variabile X
    sd_x_multiplier = 2.0,      # Dev. std. di X = sd_x_multiplier * sd(errori)
    num_estimation_iter = 50,   # Numero di iterazioni per la stima BML
    psi_initial_guess = 0.1,    # Stima iniziale per psi_BML
    verbose = FALSE             # Se TRUE, stampa i progressi dell'iterazione
) {
  
  n_total <- N_rows * N_cols # Numero totale di osservazioni
  
  # 1. Genera coordinate
  coords <- expand.grid(x = 1:N_cols, y = 1:N_rows)
  
  # 2. Genera la variabile X spazialmente correlata
  dist_matrix_x <- calculate_distance_matrix(coords)
  R_x <- create_correlation_matrix(dist_matrix_x, phi_decay_x)
  diag(R_x) <- diag(R_x) + 1e-9 # Per stabilità numerica
  
  var_x_true <- (sd_x_multiplier * sqrt(sigma_sq_true_errors))^2
  Sigma_X <- var_x_true * R_x
  diag(Sigma_X) <- diag(Sigma_X) + 1e-9 
  X_correlated <- MASS::mvrnorm(n = 1, mu = rep(0, n_total), Sigma = Sigma_X)
  X_spatial <- as.vector(X_correlated)
  
  # 3. Genera errori spazialmente correlati (epsilon)
  dist_matrix_errors <- calculate_distance_matrix(coords)
  R_errors <- create_correlation_matrix(dist_matrix_errors, phi_decay_errors)
  diag(R_errors) <- diag(R_errors) + 1e-9
  
  Sigma_errors <- sigma_sq_true_errors * R_errors
  diag(Sigma_errors) <- diag(Sigma_errors) + 1e-9
  
  errors <- MASS::mvrnorm(n = 1, mu = rep(0, n_total), Sigma = Sigma_errors)
  errors_spatial <- as.vector(errors)
  
  # 4. Genera la variabile dipendente Y
  Y_spatial <- beta_true * X_spatial + errors_spatial
  
  # 5. Centra le variabili X e Y
  X_centered <- X_spatial - mean(X_spatial)
  Y_centered <- Y_spatial - mean(Y_spatial)
  
  # 6. Codifica bivariata: seleziona coppie di vicini
  # Strategia: coppie orizzontali ((r,c), (r,c+1)) per righe dispari e colonne di partenza dispari
  paired_data <- list()
  pair_idx <- 0
  for (r in 1:N_rows) {
    if (r %% 2 == 1) { # Righe dispari
      for (c in 1:N_cols) {
        if (c %% 2 == 1 && (c + 1) <= N_cols) { # Colonne dispari, e assicura che la coppia esista
          pair_idx <- pair_idx + 1
          
          idx1 <- (r - 1) * N_cols + c       
          idx2 <- (r - 1) * N_cols + (c + 1) 
          
          paired_data[[pair_idx]] <- list(
            y_i = Y_centered[idx1], x_i = X_centered[idx1],
            y_l = Y_centered[idx2], x_l = X_centered[idx2]
          )
        }
      }
    }
  }
  
  q_pairs <- length(paired_data) 
  if (q_pairs == 0) {
    warning("Nessuna coppia formata. Controllare le dimensioni della griglia e la strategia di accoppiamento.")
    return(list(beta_hat = NA, sigma_sq_hat = NA, psi_hat = NA, q_pairs = 0))
  }
  
  # 7. Calcola le statistiche sufficienti (alpha_1 a alpha_6)
  alpha_1 <- 0; alpha_2 <- 0; alpha_3 <- 0; alpha_4 <- 0; alpha_5 <- 0; alpha_6 <- 0
  
  for (k in 1:q_pairs) {
    pair <- paired_data[[k]]
    y_i <- pair$y_i; x_i <- pair$x_i
    y_l <- pair$y_l; x_l <- pair$x_l
    
    alpha_1 <- alpha_1 + (x_i^2 + x_l^2)
    alpha_2 <- alpha_2 + (y_i^2 + y_l^2)
    alpha_3 <- alpha_3 + (x_i * y_i + x_l * y_l) 
    alpha_4 <- alpha_4 + (x_i * y_l + x_l * y_i) 
    alpha_5 <- alpha_5 + (x_i * x_l)             
    alpha_6 <- alpha_6 + (y_i * y_l)             
  }
  
  # 8. Stima BML iterativa per beta, sigma_sq, psi
  psi_hat <- psi_initial_guess
  beta_hat <- NA
  sigma_sq_hat <- NA
  
  for (iter in 1:num_estimation_iter) {
    # Stima beta_hat (Eq. 7)
    denominator_beta <- alpha_1 - 2 * psi_hat * alpha_5
    if (abs(denominator_beta) < 1e-9) {
      beta_hat <- ifelse(is.na(beta_hat), 0, beta_hat) 
    } else {
      beta_hat <- (alpha_3 - psi_hat * alpha_4) / denominator_beta
    }
    
    # Stima sigma_sq_hat (Eq. 8)
    numerator_sigma_sq <- alpha_2 + beta_hat^2 * alpha_1 - 2 * beta_hat * alpha_3 - 
      2 * psi_hat * alpha_6 - 2 * psi_hat * beta_hat^2 * alpha_5 + 
      2 * psi_hat * beta_hat * alpha_4
    denominator_sigma_sq <- 2 * q_pairs * (1 - psi_hat^2)
    
    if (abs(denominator_sigma_sq) < 1e-9 || numerator_sigma_sq < 0) {
      sigma_sq_hat <- ifelse(is.na(sigma_sq_hat) || sigma_sq_hat <=0, 1e-6, sigma_sq_hat) 
    } else {
      sigma_sq_hat <- numerator_sigma_sq / denominator_sigma_sq
      if (sigma_sq_hat <= 0) { 
        sigma_sq_hat <- 1e-6 
      }
    }
    
    # Stima psi_hat (Eq. 9)
    denominator_psi <- q_pairs * sigma_sq_hat
    if (abs(denominator_psi) < 1e-9) {
      psi_hat_new <- ifelse(iter==1 && psi_hat == psi_initial_guess, 0, psi_hat) 
    } else {
      psi_hat_new <- (alpha_6 - beta_hat * alpha_4 + beta_hat^2 * alpha_5) / denominator_psi
    }
    
    psi_hat_new <- max(-0.999, min(0.999, psi_hat_new)) # Limita psi_hat
    
    if (verbose && iter %% 10 == 0) {
      cat(paste("Iter:", iter, "beta_hat:", round(beta_hat,4), 
                "sigma_sq_hat:", round(sigma_sq_hat,4), "psi_hat:", round(psi_hat,4), "\n"))
    }
    psi_hat <- psi_hat_new
  }
  
  return(list(
    beta_hat = beta_hat,
    sigma_sq_hat = sigma_sq_hat,
    psi_hat = psi_hat,
    q_pairs = q_pairs
  ))
}

#-------------------------------------------------------------------------------
# Funzione Principale per la Simulazione Monte Carlo
#-------------------------------------------------------------------------------
monte_carlo_simulation <- function(
    num_replications = 100,     # Numero di replicazioni Monte Carlo
    N_rows = 10, N_cols = 10,   # Dimensioni della griglia
    beta_true = 1.0,
    sigma_sq_true_errors = 1.0,
    phi_decay_errors = 1.0,     # Per psi_true ~ 0.367 per d=1
    phi_decay_x = 1.0,
    sd_x_multiplier = 2.0,
    num_estimation_iter = 50,
    psi_initial_guess = 0.2
) {
  
  cat(paste("Avvio della simulazione Monte Carlo con", num_replications, "replicazioni...\n"))
  
  # Inizializza i vettori per memorizzare i risultati
  beta_estimates <- numeric(num_replications)
  sigma_sq_estimates <- numeric(num_replications)
  psi_estimates <- numeric(num_replications)
  
  # Calcola il vero psi per celle adiacenti (d=1) basato su phi_decay_errors
  psi_true_for_comparison <- exp(-(1/phi_decay_errors)^2)
  
  # Loop Monte Carlo
  for (i in 1:num_replications) {
    if (i %% (num_replications/10) == 0 || i == 1 || i == num_replications) { # Stampa il progresso
      cat(paste("Replicazione Monte Carlo:", i, "/", num_replications, "\n"))
    }
    
    # Esegui una singola simulazione e stima
    # Imposta un seme diverso per ogni replicazione per garantire la variabilità
    set.seed(Sys.time() + i) # Un modo semplice per variare il seme
    
    single_run_results <- run_spatial_simulation_and_estimation(
      N_rows = N_rows, N_cols = N_cols,
      beta_true = beta_true,
      sigma_sq_true_errors = sigma_sq_true_errors,
      phi_decay_errors = phi_decay_errors,
      phi_decay_x = phi_decay_x,
      sd_x_multiplier = sd_x_multiplier,
      num_estimation_iter = num_estimation_iter,
      psi_initial_guess = psi_initial_guess,
      verbose = FALSE # Disabilita output verboso per le singole esecuzioni durante MC
    )
    
    # Memorizza le stime
    beta_estimates[i] <- single_run_results$beta_hat
    sigma_sq_estimates[i] <- single_run_results$sigma_sq_hat
    psi_estimates[i] <- single_run_results$psi_hat
  }
  
  cat("Simulazione Monte Carlo completata.\n\n")
  
  # Rimuovi eventuali NA risultanti da problemi di stima (sebbene il codice cerchi di evitarli)
  beta_estimates <- na.omit(beta_estimates)
  sigma_sq_estimates <- na.omit(sigma_sq_estimates)
  psi_estimates <- na.omit(psi_estimates)
  
  # Calcola le statistiche riassuntive
  summary_stats <- data.frame(
    Parametro = c("beta", "sigma_sq", "psi"),
    Valore_Vero = c(beta_true, sigma_sq_true_errors, psi_true_for_comparison),
    Media_Stime = c(mean(beta_estimates), mean(sigma_sq_estimates), mean(psi_estimates)),
    Bias = c(mean(beta_estimates) - beta_true, 
             mean(sigma_sq_estimates) - sigma_sq_true_errors, 
             mean(psi_estimates) - psi_true_for_comparison),
    Bias_Relativo = c((mean(beta_estimates) - beta_true) / beta_true,
                      (mean(sigma_sq_estimates) - sigma_sq_true_errors) / sigma_sq_true_errors,
                      (mean(psi_estimates) - psi_true_for_comparison) / psi_true_for_comparison),
    MSE = c(mean((beta_estimates - beta_true)^2),
            mean((sigma_sq_estimates - sigma_sq_true_errors)^2),
            mean((psi_estimates - psi_true_for_comparison)^2))
  )
  
  # Arrotonda i risultati per una migliore visualizzazione
  summary_stats[, 3:6] <- round(summary_stats[, 3:6], 5)
  
  print(summary_stats)
  
  return(list(
    summary_statistics = summary_stats,
    beta_estimates = beta_estimates,
    sigma_sq_estimates = sigma_sq_estimates,
    psi_estimates = psi_estimates
  ))
}

#-------------------------------------------------------------------------------
# Esempio di Utilizzo della Simulazione Monte Carlo
#-------------------------------------------------------------------------------
# Imposta il seme principale per la riproducibilità della simulazione MC nel suo complesso
# (le singole esecuzioni all'interno usano semi variabili)
set.seed(456) 

# Esegui la simulazione Monte Carlo
# Questi parametri sono simili all'Esperimento 1 nella Tabella 1 del paper
# (n=100, psi=0.367).
# phi_decay_errors = 1.0 dà psi_true ~ 0.367 per d=1.
# sigma_sq_true_errors = 1.0 è il valore vero della varianza degli errori nel paper.
monte_carlo_results <- monte_carlo_simulation(
  num_replications = 1000,      # Ridotto per esecuzione rapida, aumentare per risultati migliori
  N_rows = 10, N_cols = 10,   # n = 100
  beta_true = 1.0,
  sigma_sq_true_errors = 1.0, 
  phi_decay_errors = 1.0,     
  phi_decay_x = 1.0,          # Come nel paper
  sd_x_multiplier = 2.0,      # Come nel paper
  num_estimation_iter = 50,
  psi_initial_guess = 0
)

# Per visualizzare le distribuzioni delle stime (opzionale):
# par(mfrow=c(1,3))
# hist(monte_carlo_results$beta_estimates, main="Distribuzione di beta_hat", xlab="beta_hat")
# abline(v = 1.0, col="red", lwd=2)
# hist(monte_carlo_results$sigma_sq_estimates, main="Distribuzione di sigma_sq_hat", xlab="sigma_sq_hat")
# abline(v = 1.0, col="red", lwd=2)
# hist(monte_carlo_results$psi_estimates, main="Distribuzione di psi_hat", xlab="psi_hat")
# abline(v = exp(-(1/1.0)^2), col="red", lwd=2) # psi_true per phi_decay_errors=1.0
# par(mfrow=c(1,1))


#-----------------------------------------------------------------------------#
# 0. Caricamento Librerie e Impostazioni Iniziali
#-----------------------------------------------------------------------------#
library(sf)
library(spdep)
library(spatialreg) # Per GMerrorsar
library(Matrix)     # Per matrici sparse CsparseMatrix
library(ggplot2)    # Opzionale, per visualizzare i dati o i residui
library(dplyr)      # Per la manipolazione dei dati

set.seed(12345) # Per riproducibilità

#-----------------------------------------------------------------------------#
# 1. Funzione per Generare u con Approssimazione Taylor/Neumann
#-----------------------------------------------------------------------------#
generate_u_taylor_approx <- function(epsilon_vec, W_sparse_matrix, lambda, K_max, verbose = TRUE) {
  if (abs(lambda) < 1e-9 || K_max == 0) { # Se lambda è zero o non ci sono termini, u = epsilon
    if (verbose) cat("  Approssimazione Taylor: lambda=0 o K_max=0, u = epsilon.\n")
    return(epsilon_vec)
  }
  
  N <- length(epsilon_vec)
  u_approx <- epsilon_vec # Termine per k=0: (lambda*W)^0 * epsilon = epsilon
  
  current_power_of_lambdaW_times_epsilon <- epsilon_vec # Questo è ((lambda*W)^0)*epsilon
  
  if (verbose) cat("  Approssimazione Taylor: ")
  for (k_iter in 1:K_max) {
    if (verbose && (k_iter %% 10 == 0 || k_iter == K_max)) cat(paste0(k_iter, "..."))
    
    # Calcola (W * precedente_termine_componente_epsilon)
    # poi moltiplica per lambda per ottenere il termine (lambda*W)^k_iter * epsilon
    current_power_of_lambdaW_times_epsilon <- lambda * (W_sparse_matrix %*% current_power_of_lambdaW_times_epsilon)
    
    u_approx <- u_approx + current_power_of_lambdaW_times_epsilon
    
    # Opzionale: Criterio di stop se l'ultimo termine aggiunto è molto piccolo
    # if (sum(current_power_of_lambdaW_times_epsilon^2) < 1e-16 * sum(epsilon_vec^2)) { # Tolleranza relativa
    #   if (verbose) cat(paste0("\n  Convergenza approssimazione Taylor raggiunta a k=", k_iter, "\n"))
    #   break
    # }
  }
  if (verbose) cat(" Fatto.\n")
  return(u_approx)
}

#-----------------------------------------------------------------------------#
# 2. Funzione per Generare Dati Spaziali di Test con u Approssimato
#-----------------------------------------------------------------------------#
generate_spatial_data_approx_dgp <- function(N_obs = 10000, 
                                             beta_true_vec = c(1.0, 1.5), # Intercetta, beta1
                                             lambda_dgp = 0.7,      # Lambda "vero" per il DGP approssimato
                                             sigma_sq_eps_dgp = 1.0, # Varianza rumore bianco per il DGP
                                             k_neighbors_w = 6,      # Per la matrice W
                                             K_max_taylor = 30,     # Termini per l'approssimazione Taylor
                                             use_exact_invIrW_for_comparison = FALSE, # Per confronto se N è piccolo
                                             verbose_approx = TRUE) {
  
  cat(paste0("Inizio generazione dati (DGP approssimato) per N = ", N_obs, "...\n"))
  
  # Genera coordinate casuali
  coords <- data.frame(x_coord = runif(N_obs), y_coord = runif(N_obs))
  points_sf <- st_as_sf(coords, coords = c("x_coord", "y_coord"), crs = 4326)
  
  # Crea covariata X
  X_mat <- cbind(Intercept = beta_true_vec[1], x_value = rnorm(N_obs, mean = 0, sd = 1))
  # Modificato per usare beta_true_vec[1] come intercetta fissa nella X_mat per il DGP
  # e beta_true_vec[2] come coefficiente per x_value
  mean_signal_fixed_part <- X_mat[,1] # Solo intercetta
  mean_signal_X_part <- X_mat[,2] * beta_true_vec[2] # x_value * beta1
  
  # Crea matrice di pesi W (k-nearest neighbors)
  cat("  Creazione matrice di pesi W...\n")
  knn <- knearneigh(st_coordinates(points_sf), k = k_neighbors_w)
  nb <- knn2nb(knn)
  listw_obj <- nb2listw(nb, style = "W", zero.policy = TRUE) # W normalizzata per riga
  W_mat_sparse <- as(listw_obj, "CsparseMatrix")
  
  # Genera errori epsilon i.i.d.
  epsilon <- rnorm(N_obs, mean = 0, sd = sqrt(sigma_sq_eps_dgp))
  
  # Genera errori spazialmente correlati u
  cat(paste0("  Generazione u con approssimazione Taylor (K_max=", K_max_taylor, ")...\n"))
  time_start_u_approx <- Sys.time()
  u_approximated <- generate_u_taylor_approx(epsilon_vec = epsilon, 
                                             W_sparse_matrix = W_mat_sparse, 
                                             lambda = lambda_dgp, 
                                             K_max = K_max_taylor,
                                             verbose = verbose_approx)
  time_end_u_approx <- Sys.time()
  time_gen_u_approx <- difftime(time_end_u_approx, time_start_u_approx, units="secs")
  cat(paste0("    Tempo generazione u (approssimato): ", round(as.numeric(time_gen_u_approx),3), " sec\n"))
  
  u_exact <- NULL
  time_gen_u_exact <- NA
  if (use_exact_invIrW_for_comparison) {
    cat("  Generazione u con invIrW esatto (per confronto)...\n")
    time_start_u_exact <- Sys.time()
    u_exact <- tryCatch({
      inv_I_lambda_W <- invIrW(W_mat_sparse, rho = lambda_dgp)
      as.vector(inv_I_lambda_W %*% epsilon)
    }, error = function(e) { 
      warning(paste("Errore in invIrW:", e$message)); NULL 
    })
    time_end_u_exact <- Sys.time()
    if(!is.null(u_exact)) {
      time_gen_u_exact <- difftime(time_end_u_exact, time_start_u_exact, units="secs")
      cat(paste0("    Tempo generazione u (esatto): ", round(as.numeric(time_gen_u_exact),3), " sec\n"))
      
      # Confronto rapido (opzionale)
      if(length(u_approximated) == length(u_exact)) {
        abs_diff <- abs(u_approximated - u_exact)
        cat(paste0("    Confronto u_approx vs u_exact: MSE = ", mean(abs_diff^2), 
                   ", Max Abs Diff = ", max(abs_diff), "\n"))
      }
    }
  }
  
  # Genera variabile dipendente Y usando u_approximated
  # y_vec <- as.vector(X_mat %*% beta_true_vec + u_approximated) # Se X_mat ha Intercept e x_value
  y_vec <- as.vector(mean_signal_fixed_part + mean_signal_X_part + u_approximated)
  
  
  data_sf <- points_sf
  data_sf$y_gm <- y_vec
  data_sf$x_gm <- X_mat[, "x_value"] 
  
  cat(paste0("Dati (DGP approssimato) generati (N=", N_obs, ").\n"))
  return(list(data = data_sf, 
              listw = listw_obj, 
              params_dgp = list(beta0 = beta_true_vec[1], beta1 = beta_true_vec[2], 
                                lambda = lambda_dgp, sigma_sq_eps = sigma_sq_eps_dgp),
              u_approximated = u_approximated, 
              u_exact = u_exact, # Sarà NULL se non calcolato
              time_gen_u_approx_sec = as.numeric(time_gen_u_approx),
              time_gen_u_exact_sec = if(!is.null(u_exact)) as.numeric(time_gen_u_exact) else NA
  ))
}

#-----------------------------------------------------------------------------#
# 3. Parametri per l'Esperimento
#-----------------------------------------------------------------------------#
N_experiment <- 10000    # Numero di osservazioni 
# (ATTENZIONE: use_exact_invIrW_for_comparison = TRUE con N>5000 può essere MOLTO lento)
K_taylor_terms <- 30      # Numero di termini per l'approssimazione di Taylor
lambda_dgp_true <- 0.7
sigma_sq_eps_dgp_true <- 1.0
beta_dgp_true <- c(1.0, 1.5) # Intercetta, beta1 per x_gm

# Genera i dati una volta usando l'approssimazione
cat(paste0("--- Inizio Esperimento: Generazione Dati (N=", N_experiment, ", K_taylor=", K_taylor_terms, ") ---\n"))
dgp_results <- generate_spatial_data_approx_dgp(
  N_obs = N_experiment,
  beta_true_vec = beta_dgp_true,
  lambda_dgp = lambda_dgp_true,
  sigma_sq_eps_dgp = sigma_sq_eps_dgp_true,
  k_neighbors_w = 6,
  K_max_taylor = K_taylor_terms,
  use_exact_invIrW_for_comparison = ifelse(N_experiment <= 3000, TRUE, FALSE), # Confronta solo se N non troppo grande
  verbose_approx = TRUE
)

df_experiment <- dgp_results$data
listw_experiment <- dgp_results$listw

#-----------------------------------------------------------------------------#
# 4. Stima con GMerrorsar sui Dati Approssimati
#-----------------------------------------------------------------------------#
cat("\n--- Stima Modello GMerrorsar sui Dati Generati con Approssimazione ---\n")
formula_exp <- y_gm ~ x_gm

time_start_gmm <- Sys.time()
model_gmm_on_approx_data <- tryCatch({
  GMerrorsar(
    formula = formula_exp, 
    data = df_experiment, 
    listw = listw_experiment,
    zero.policy = TRUE, 
    verbose = FALSE 
  )
}, error = function(e) {
  cat(paste0("ERRORE durante la stima GMerrorsar: ", e$message, "\n"))
  NULL
})
time_end_gmm <- Sys.time()
time_estimation_gmm <- difftime(time_end_gmm, time_start_gmm, units="secs")

#-----------------------------------------------------------------------------#
# 5. Risultati
#-----------------------------------------------------------------------------#
cat("\n\n--- RIEPILOGO RISULTATI ---\n")
cat(paste0("Parametri del DGP (usati per generare dati con u approssimato):\n"))
cat(paste0("  Beta0 Vero: ", dgp_results$params_dgp$beta0, "\n"))
cat(paste0("  Beta1 Vero: ", dgp_results$params_dgp$beta1, "\n"))
cat(paste0("  Lambda Vero (DGP): ", dgp_results$params_dgp$lambda, "\n"))
cat(paste0("  Sigma^2 Eps Vero (DGP): ", dgp_results$params_dgp$sigma_sq_eps, "\n"))
cat(paste0("Tempo generazione u (approssimato, K_max=", K_taylor_terms, "): ", 
           round(dgp_results$time_gen_u_approx_sec, 3), " sec\n"))

if (!is.na(dgp_results$time_gen_u_exact_sec)) {
  cat(paste0("Tempo generazione u (esatto, invIrW): ", 
             round(dgp_results$time_gen_u_exact_sec, 3), " sec\n"))
}

cat("\nStime da GMerrorsar (su dati con u approssimato):\n")
if (!is.null(model_gmm_on_approx_data)) {
  summary_gmm <- summary(model_gmm_on_approx_data)
  cat(paste0("  Beta0 Stimato (GM): ", round(coef(model_gmm_on_approx_data)[1], 4), "\n")) # Assumendo intercetta sia il primo
  cat(paste0("  Beta1 Stimato (GM): ", round(coef(model_gmm_on_approx_data)["x_gm"], 4), "\n"))
  cat(paste0("  Lambda Stimato (GM): ", round(model_gmm_on_approx_data$lambda, 4), "\n"))
  cat(paste0("  Sigma^2 Stimata (GM): ", round(model_gmm_on_approx_data$s2, 4), "\n"))
  cat(paste0("  Tempo stima GMerrorsar: ", round(as.numeric(time_estimation_gmm), 3), " sec\n"))
  
  # Calcolo Bias Semplice
  bias_beta0_gmm <- coef(model_gmm_on_approx_data)[1] - dgp_results$params_dgp$beta0
  bias_beta1_gmm <- coef(model_gmm_on_approx_data)["x_gm"] - dgp_results$params_dgp$beta1
  bias_lambda_gmm <- model_gmm_on_approx_data$lambda - dgp_results$params_dgp$lambda
  bias_sigma_gmm <- model_gmm_on_approx_data$s2 - dgp_results$params_dgp$sigma_sq_eps
  
  cat("\nBias delle stime GM (Stima - Vero_DGP_Approssimato):\n")
  cat(paste0("  Bias Beta0: ", round(bias_beta0_gmm, 4), "\n"))
  cat(paste0("  Bias Beta1: ", round(bias_beta1_gmm, 4), "\n"))
  cat(paste0("  Bias Lambda: ", round(bias_lambda_gmm, 4), "\n"))
  cat(paste0("  Bias Sigma^2: ", round(bias_sigma_gmm, 4), "\n"))
  
} else {
  cat("  Stima GMerrorsar fallita.\n")
}

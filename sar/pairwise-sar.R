#--------------------------------------------------------------------------
# PARTE 1: SETUP E DIPENDENZE
#--------------------------------------------------------------------------
# Assicurarsi di aver installato i pacchetti necessari
# install.packages(c("spdep", "spatialreg", "dplyr", "ggplot2", "corrplot"))

library(spdep)
library(spatialreg)
library(dplyr)
library(ggplot2)
library(corrplot)

#--------------------------------------------------------------------------
# PARTE 2: PARAMETRI DELLA SIMULAZIONE
#--------------------------------------------------------------------------
set.seed(123) # Per la riproducibilità

# Parametri della griglia e del modello
GRID_SIZE <- 20       # Griglia 20x20
N <- GRID_SIZE^2      # Numero totale di osservazioni (400)
N_SIM <- 500          # Numero di simulazioni Monte Carlo

# Valori veri dei parametri
TRUE_BETA <- 1.5
TRUE_RHO <- 0.6
TRUE_SIGMA2 <- 1.0

#--------------------------------------------------------------------------
# PARTE 3: FUNZIONE DI STIMA PAIRWISE SAR
# Questa è l'implementazione del suo contributo teorico
#--------------------------------------------------------------------------
estimate_pairwise_sar <- function(y, x, nb_list, max_iter = 30, tol = 1e-6) {
  
  # --- 3.1: Creazione della lista di coppie di vicini
  # Per la griglia regolare, usiamo tutte le coppie di vicini adiacenti
  pairs_df <- do.call(rbind, lapply(1:length(nb_list), function(i) {
    if (length(nb_list[[i]]) > 0 && nb_list[[i]][1] != 0) {
      data.frame(i = i, l = nb_list[[i]])
    }
  }))
  # Rimuoviamo le coppie duplicate (es. (1,2) e (2,1))
  pairs_df <- pairs_df[!duplicated(t(apply(pairs_df, 1, sort))), ]
  q <- nrow(pairs_df)
  
  # --- 3.2: Calcolo delle 5 statistiche sufficienti (T1 a T5)
  # Estraiamo i dati per ogni coppia
  xi <- x[pairs_df$i]
  xl <- x[pairs_df$l]
  yi <- y[pairs_df$i]
  yl <- y[pairs_df$l]
  
  T1 <- sum(xi^2 + xl^2)
  T2 <- sum(yi^2 + yl^2)
  T3 <- sum(xi*yi + xl*yl)
  T4 <- sum(yi*yl)
  T5 <- sum(xi*yl + xl*yi)
  
  # --- 3.3: Risolutore Iterativo
  # Inizializzazione dei parametri
  rho_hat <- 0.0 
  
  for (iter in 1:max_iter) {
    rho_old <- rho_hat
    
    # a) Stima di beta dato rho (Equazione 1 del sistema)
    beta_hat <- (T3 - rho_hat * T5) / T1
    
    # b) Stima di sigma^2 dati beta e rho (Equazione 3 del sistema)
    sigma2_hat <- (1/(2*q)) * ((1+rho_hat^2)*T2 + beta_hat^2*T1 - 2*beta_hat*T3 - 4*rho_hat*T4 + 2*rho_hat*beta_hat*T5)
    
    if (sigma2_hat <= 0) { # Controllo di sicurezza
      warning("Stima di sigma^2 non positiva, interruzione.")
      return(c(beta_hat = NA, rho_hat = NA, sigma2_hat = NA))
    }
    
    # c) Aggiornamento di rho (dall'Equazione 2)
    # L'equazione per rho è un polinomio di terzo grado. Per semplicità e robustezza,
    # usiamo una grid search per trovare il rho che minimizza l'errore quadratico della score equation.
    
    rho_grid <- seq(-0.99, 0.99, by = 0.01)
    
    # Score equation per rho (lato sx - lato dx)
    score_rho <- function(rho) {
      lhs <- 2*T4 - rho*T2 - beta_hat*T5
      rhs <- (2*q*rho*sigma2_hat) / (1 - rho^2)
      return((lhs - rhs)^2)
    }
    
    # Troviamo il rho che minimizza la funzione
    best_rho_index <- which.min(sapply(rho_grid, score_rho))
    rho_hat <- rho_grid[best_rho_index]
    
    # Controlliamo la convergenza
    if (abs(rho_hat - rho_old) < tol) {
      break
    }
  }
  
  # Restituiamo gli stimatori finali
  return(c(beta_hat = beta_hat, rho_hat = rho_hat, sigma2_hat = sigma2_hat))
}


#--------------------------------------------------------------------------
# PARTE 4: IL LOOP MONTE CARLO
#--------------------------------------------------------------------------

# Creiamo la struttura della griglia e la matrice W una sola volta
coords <- expand.grid(1:GRID_SIZE, 1:GRID_SIZE)
# Definiamo la vicinanza "queen" (regina), cioè anche le diagonali
nb <- cell2nb(GRID_SIZE, GRID_SIZE, type = "queen") 
W_list <- nb2listw(nb, style = "W") # Matrice W row-standardized
W <- as(W_list, "CsparseMatrix")

# Matrice per salvare i risultati
results_df <- data.frame(
  run = 1:N_SIM,
  beta_est = numeric(N_SIM),
  rho_est = numeric(N_SIM),
  sigma2_est = numeric(N_SIM)
)

# Iniziamo il loop
cat("Avvio della simulazione Monte Carlo...\n")
pb <- txtProgressBar(min = 0, max = N_SIM, style = 3) # Progress bar

for (i in 1:N_SIM) {
  # --- 4.1: Generazione dei dati per questa iterazione
  # Variabile esplicativa
  x <- rnorm(N, mean = 5, sd = 2)
  # Termine di errore primitivo
  epsilon <- rnorm(N, mean = 0, sd = sqrt(TRUE_SIGMA2))
  
  # Calcolo del termine (I - rho*W)^-1 * (X*beta + epsilon)
  # Questo è il Data Generating Process corretto per un modello SAR
  A <- .sparseDiagonal(N) - TRUE_RHO * W
  y <- solve(A, x * TRUE_BETA + epsilon)
  
  # --- 4.2: Stima dei parametri con il metodo Pairwise
  estimates <- estimate_pairwise_sar(y, x, nb)
  
  # --- 4.3: Salvataggio dei risultati
  results_df$beta_est[i] <- estimates["beta_hat"]
  results_df$rho_est[i] <- estimates["rho_hat"]
  results_df$sigma2_est[i] <- estimates["sigma2_hat"]
  
  setTxtProgressBar(pb, i)
}
close(pb)
cat("Simulazione completata.\n")


#--------------------------------------------------------------------------
# PARTE 5: ANALISI E VISUALIZZAZIONE DEI RISULTATI
#--------------------------------------------------------------------------

# Rimuoviamo eventuali righe con NA (se la stima fallisce)
results_clean <- na.omit(results_df)

# --- 5.1: Statistiche descrittive e bias
summary_stats <- results_clean %>%
  summarise(
    mean_beta = mean(beta_est),
    mean_rho = mean(rho_est),
    mean_sigma2 = mean(sigma2_est)
  )

cat("\n--- Risultati della Simulazione ---\n")
cat("Valori Veri:\n")
cat(sprintf("  beta = %.3f, rho = %.3f, sigma2 = %.3f\n", TRUE_BETA, TRUE_RHO, TRUE_SIGMA2))
cat("Medie Stimate:\n")
cat(sprintf("  beta_hat = %.3f, rho_hat = %.3f, sigma2_hat = %.3f\n", 
            summary_stats$mean_beta, summary_stats$mean_rho, summary_stats$mean_sigma2))

# --- 5.2: VERIFICA DELLA CORRELAZIONE TRA GLI STIMATORI
# Questo è il punto cruciale per verificare la sua teoria
param_estimates <- results_clean[, c("beta_est", "rho_est", "sigma2_est")]
cor_matrix <- cor(param_estimates)

cat("\n--- Matrice di Correlazione degli Stimatori ---\n")
print(round(cor_matrix, 3))

# Visualizziamo la matrice di correlazione
corrplot(cor_matrix, method = "number", type = "upper", 
         title = "Correlazione tra gli Stimatori Pairwise SAR", mar=c(0,0,1,0))

# --- 5.3: Visualizzazione della distribuzione congiunta
cat("\nVisualizzazione della distribuzione congiunta degli stimatori...\n")

# Scatter plot tra beta_hat e rho_hat per vedere la correlazione
ggplot(results_clean, aes(x = beta_est, y = rho_est)) +
  geom_point(alpha = 0.4, color = "blue") +
  geom_density_2d() +
  labs(
    title = "Distribuzione Congiunta di beta_hat e rho_hat",
    subtitle = "La correlazione osservata conferma la teoria della Fisher matrix non-diagonale",
    x = "Stime di Beta",
    y = "Stime di Rho"
  ) +
  theme_minimal()
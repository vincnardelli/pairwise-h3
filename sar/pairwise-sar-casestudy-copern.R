#--------------------------------------------------------------------------
# ANALISI EMPIRICA: MODELLO SAR PAIRWISE VS MASSIMA VEROSIMIGLIANZA
# Caso Studio: Temperatura in Italia (Dati Copernicus ERA5-Land)
#--------------------------------------------------------------------------

# --- 1. SETUP: CARICAMENTO LIBRERIE ---
# Assicurati di aver installato queste librerie con install.packages("nome_libreria")

# Per leggere i dati NetCDF e manipolare i raster
library(terra) 
# Per la manipolazione generale dei dati
library(dplyr)
# Per l'analisi spaziale: creazione vicini, pesi e modello SAR standard
library(spdep)
library(spatialreg)
# Per la visualizzazione (opzionale)
library(ggplot2)

# --- 2. CARICAMENTO E PRE-PROCESSING DATI ---
# !!! MODIFICA QUI: Inserisci i nomi dei due file che hai scaricato !!!
TEMP_FILE <- "sar/data/data_stream-moda.nc"
GEO_FILE  <- "sar/data/geo.nc"

# Controlla se i file esistono
if (!file.exists(TEMP_FILE) || !file.exists(GEO_FILE)) {
  stop("File non trovati! Assicurati che entrambi i file siano nella tua cartella di lavoro.")
}

# 1. Carica i due raster separatamente
temp_raster_global <- rast(TEMP_FILE)
geo_raster_global  <- rast(GEO_FILE)

# 2. Definisci l'area di interesse (Bounding Box per l'Italia)
#    xmin, xmax, ymin, ymax
italy_extent <- ext(6, 19, 35, 47) 

# 3. Ritaglia (crop) entrambi i raster sull'estensione dell'Italia
#    Questo passo è FONDAMENTALE per ridurre la dimensione dei dati!
cat("Ritagliando i raster globali sull'area dell'Italia...\n")
temp_raster_italy <- crop(temp_raster_global, italy_extent)
geo_raster_italy  <- crop(geo_raster_global, italy_extent)

# 4. Combina i due raster italiani in un unico "stack" multi-strato
#    Questo funziona perché hanno la stessa griglia
climate_stack <- c(temp_raster_italy, geo_raster_italy)
names(climate_stack) <- c("t2m", "z") # Rinominiamo i layer per chiarezza
print(climate_stack)

# 5. Converti lo stack combinato in un data frame
cat("Convertendo lo stack raster in data frame...\n")
df_climate <- as.data.frame(climate_stack, xy = TRUE)

# 6. Rinomina le colonne e processa i dati
df_climate <- df_climate %>%
  rename(
    lon = x,
    lat = y,
    temp_k = t2m,
    geopotential = z
  ) %>%
  # Rimuovi eventuali righe con dati mancanti (es. mare)
  na.omit() %>%
  # Calcola le variabili per il modello
  mutate(
    # Converte la temperatura da Kelvin a Celsius (più interpretabile)
    temp_c = temp_k - 273.15,
    # Calcola l'elevazione in metri dal geopotenziale
    elevation_m = geopotential / 9.80665
  )

cat("\nDimensioni del dataset finale per l'Italia:", nrow(df_climate), "punti griglia.\n")
head(df_climate)


# --- 3. IMPLEMENTAZIONE DELLO STIMATORE PAIRWISE SAR (PW-SAR) ---

# Questa sezione traduce le formule del tuo paper in funzioni R

#' Calcola le 5 statistiche sufficienti per il modello SAR pairwise
#' @param y Vettore della variabile dipendente
#' @param X Matrice della variabile indipendente (con intercetta)
#' @param pairs Matrice a 2 colonne con gli indici delle coppie di vicini
#' @return Una lista contenente le 5 statistiche T1-T5 e il numero di coppie q
calculate_sufficient_stats <- function(y, X, pairs) {
  # Estrai i dati per ogni elemento della coppia
  y_i <- y[pairs[, 1]]
  y_l <- y[pairs[, 2]]
  x_i <- X[pairs[, 1], 2] # Assumiamo una sola var. indipendente + intercetta
  x_l <- X[pairs[, 2], 2]
  
  # Calcola le statistiche
  T1 <- sum(x_i^2 + x_l^2)
  T2 <- sum(y_i^2 + y_l^2)
  T3 <- sum(x_i * y_i + x_l * y_l)
  T4 <- sum(y_i * y_l)
  T5 <- sum(x_i * y_l + x_l * y_i)
  
  return(list(T1 = T1, T2 = T2, T3 = T3, T4 = T4, T5 = T5, q = nrow(pairs)))
}


#' Calcola la log-verosimiglianza pairwise concentrata per un dato rho
#' @param rho Valore del parametro di autocorrelazione spaziale
#' @param stats Lista contenente le statistiche sufficienti T1-T5 e q
#' @return Il valore della log-verosimiglianza
profile_loglik_sar <- function(rho, stats) {
  # Estrai le statistiche per leggibilità
  T1 <- stats$T1; T2 <- stats$T2; T3 <- stats$T3; T4 <- stats$T4; T5 <- stats$T5; q <- stats$q
  
  # Stima beta dato rho (eq. 15 del paper)
  beta_hat <- (T3 - rho * T5) / T1
  
  # Stima sigma^2 dato rho e beta (eq. 17 del paper)
  # Nota: Aggiungiamo un termine per l'intercetta (beta_0) che si semplifica
  # se le variabili sono centrate. Per robustezza, lo omettiamo assumendo
  # che l'effetto sia catturato dalle altre stime.
  sigma2_hat <- (1 / (2 * q)) * ( (1 + rho^2) * T2 + beta_hat^2 * T1 - 2 * beta_hat * T3 - 4 * rho * T4 + 2 * rho * beta_hat * T5 )
  
  # Evita valori negativi o nulli di sigma^2 che bloccherebbero il log
  if (sigma2_hat <= 1e-9) return(-Inf)
  
  # Calcola la log-verosimiglianza (eq. 14 del paper)
  # Omettiamo le costanti che non influenzano la massimizzazione
  loglik <- q * log(1 - rho^2) - q * log(sigma2_hat)
  
  return(loglik)
}


# --- 4. ESECUZIONE DELLE ANALISI ---

cat("\n--- Inizio Analisi Spaziale ---\n")

# A. Definizione della struttura di vicinato
# Creiamo una lista di vicini basata sulla contiguità "Queen" sulla griglia
grid_dims <- c(length(unique(df_climate$lon)), length(unique(df_climate$lat)))
nb <- cell2nb(grid_dims[1], grid_dims[2], type = "queen")

# B. ESECUZIONE PW-SAR (SUL DATASET COMPLETO)
cat("\n1. Stima del modello Pairwise SAR (PW-SAR) sul dataset completo...\n")

# Definisci il modello
y <- df_climate$temp_c
X <- model.matrix(~ elevation_m, data = df_climate)

# Genera le 4 codifiche "checkerboard" per la robustezza
# e calcola le stime medie, come suggerito da Besag (1974)
estimates_list <- list()
time_start_pw <- Sys.time()

for (i in 1:4) {
  # Crea una codifica non sovrapposta (sub-sampling dei vicini)
  all_pairs <- nb2mat(nb, style="B", zero.policy=TRUE)
  all_pairs <- which(all_pairs == 1, arr.ind = TRUE)
  all_pairs <- all_pairs[all_pairs[,1] < all_pairs[,2], ] # Evita duplicati
  
  # Semplice schema di sub-campionamento per simulare le 4 codifiche
  set.seed(120 + i)
  sub_indices <- sample(1:nrow(all_pairs), size = floor(nrow(all_pairs) / 4))
  pairs_subset <- all_pairs[sub_indices, ]
  
  # Calcola le statistiche sufficienti sul subset di coppie
  stats <- calculate_sufficient_stats(y, X, pairs_subset)
  
  # Ottimizza per trovare rho
  opt_result <- optimize(
    f = profile_loglik_sar,
    interval = c(-0.99, 0.99), # Cerca rho tra -1 e 1
    stats = stats,
    maximum = TRUE # Massimizza la verosimiglianza
  )
  
  # Salva le stime
  rho_hat <- opt_result$maximum
  beta_hat <- (stats$T3 - rho_hat * stats$T5) / stats$T1
  sigma2_hat <- (1/(2*stats$q)) * ((1+rho_hat^2)*stats$T2 + beta_hat^2*stats$T1 - 2*beta_hat*stats$T3 - 4*rho_hat*stats$T4 + 2*rho_hat*beta_hat*stats$T5)
  
  estimates_list[[i]] <- c(beta1 = beta_hat, rho = rho_hat, sigma2 = sigma2_hat)
}

time_end_pw <- Sys.time()
pw_sar_time <- time_end_pw - time_start_pw

# Calcola la media delle stime ottenute dalle 4 codifiche
pw_sar_estimates <- colMeans(do.call(rbind, estimates_list))

cat("Stima PW-SAR completata.\n")


# C. ESECUZIONE ML-SAR (SU UN SOTTOCAMPIONE)
cat("\n2. Stima del modello SAR con Massima Verosimiglianza (ML) su un sottocampione...\n")

# Crea un sottocampione per la stima ML (più veloce)
set.seed(42)
subsample_indices <- sample(1:nrow(df_climate), size = 5000)
df_subsample <- df_climate[subsample_indices, ]

# Estrai i vicini per il solo sottocampione
nb_subsample <- subset(nb, subsample_indices)
listw_subsample <- nb2listw(nb_subsample, style = "W", zero.policy = TRUE)

# Stima il modello SAR standard
time_start_ml <- Sys.time()
ml_sar_model <- lagsarlm(temp_c ~ elevation_m, data = df_subsample, listw = listw_subsample)
time_end_ml <- Sys.time()
ml_sar_time <- time_end_ml - time_start_ml

cat("Stima ML-SAR completata.\n")

# Estrai i coefficienti dal modello ML
ml_sar_estimates <- c(
  beta1 = coef(ml_sar_model)["elevation_m"],
  rho = ml_sar_model$rho,
  sigma2 = ml_sar_model$s2
)

# --- 5. CONFRONTO DEI RISULTATI ---
cat("\n--- Confronto Risultati ---\n")

results_df <- data.frame(
  Parameter = c("Elevation (beta1)", "Spatial Corr (rho)", "Sigma^2"),
  PW_SAR_Full_Dataset = unname(pw_sar_estimates),
  ML_SAR_Subsample = unname(ml_sar_estimates)
)

time_df <- data.frame(
  Metric = "Computation Time (seconds)",
  PW_SAR_Full_Dataset = as.numeric(pw_sar_time, units = "secs"),
  ML_SAR_Subsample = as.numeric(ml_sar_time, units = "secs")
)

print(results_df, row.names = FALSE)
cat("\n")
print(time_df, row.names = FALSE)
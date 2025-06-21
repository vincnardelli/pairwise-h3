# Pacchetti necessari
# install.packages(c("spatialreg", "spdep", "sf")) # Eseguire se non già installati
library(spatialreg)
library(spdep)
library(sf)

# --- 1. Preparazione dei dati e della matrice dei pesi ---
# Dati di esempio: Columbus, OH crime data (da spdep)
data(columbus, package = "spdep") # Carica il data.frame 'columbus'

# Converti il data.frame in un oggetto sf per creare i vicini
# Usiamo le coordinate X, Y presenti nel dataset e un CRS di esempio (WGS84)
# Assicurati che i nomi delle colonne delle coordinate siano corretti per il tuo dataset 'columbus'
# Tipicamente sono "X" e "Y" o "COL.extract.X", "COL.extract.Y"
# Per la versione standard di 'columbus' in 'spdep', sono 'X' e 'Y'.
col_sf <- st_as_sf(columbus, coords = c("X", "Y"), crs = 4326)

# Creazione della matrice dei pesi spaziali (es. Queen contiguity)
# Usiamo i rownames originali del dataframe 'columbus' per coerenza degli ID
col_nb <- poly2nb(col_sf, row.names = rownames(columbus))
listw_col <- nb2listw(col_nb, style = "W", zero.policy = TRUE) # zero.policy=TRUE per gestire isole

# --- 2. Definizione dei parametri della "simulazione" ---
# Valore "obiettivo" per lambda (ipotetico, per il calcolo del bias)
# In una vera simulazione Monte Carlo, questo sarebbe il valore vero usato per generare i dati.
valore_obiettivo_lambda <- 0.4 # Esempio

# Numero di "simulazioni" (qui, fit ripetuti sugli stessi dati)
n_simulazioni <- 5 # Poche run sono sufficienti per dimostrazione se i dati non cambiano

mse_lambda_runs <- numeric(n_simulazioni) # Vettore per conservare gli MSE di ogni run

set.seed(123) # Per riproducibilità

cat("Inizio 'simulazioni' (fit ripetuti su stessi dati)...\n")

# Definiamo la formula del modello usando i nomi delle colonne di 'columbus'
# Ad esempio, modelliamo CRIME in funzione di INC (income)
formula_model <- CRIME ~ INC

# --- 3. Loop di "Simulazione" ---
for (i in 1:n_simulazioni) {
  # Fit del modello GMerrorsar
  # Useremo sempre gli stessi dati (il data.frame 'columbus') e la stessa listw
  # Questo significa che stima_lambda e se_lambda saranno probabilmente identici
  # in ogni "simulazione", a meno di non convergenze stocastiche (improbabile con method="nlminb").
  
  gm_model_fit <- NULL # Inizializza a NULL
  tryCatch({
    # Il modello specificato dall'utente era:
    # GMerrorsar(y_gm~x_gm,data=df_for_gm_sf_filtered,listw=listw_gm,zero.policy=TRUE,method="nlminb")
    # Lo adattiamo al nostro dataset:
    gm_model_fit <- GMerrorsar(
      formula = formula_model,
      data = columbus, # GMerrorsar si aspetta un data.frame
      listw = listw_col,
      zero.policy = TRUE, # Come nell'esempio utente
      method = "nlminb"    # Come nell'esempio utente (è il default)
    )
  }, error = function(e) {
    cat("Errore nel fit del modello GMerrorsar per la run ", i, ": ", e$message, "\n")
    # gm_model_fit rimarrà NULL
  })
  
  if (!is.null(gm_model_fit)) {
    summary_gm <- summary(gm_model_fit) # L'oggetto summary è di classe "summary.sarlm"
    
    # --- 4. Estrazione stima di lambda e suo errore standard ---
    stima_lambda <- summary_gm$lambda   # Stima del parametro lambda
    se_lambda <- summary_gm$lambda.se # Errore standard di lambda
    
    if (is.null(stima_lambda) || is.null(se_lambda)) {
      cat("Lambda o SE di Lambda non trovati per la run ", i, ". Controllare l'output del summary.\n")
      # print(str(summary_gm)) # Decommentare per ispezionare la struttura di summary_gm
      mse_lambda_runs[i] <- NA # Segna come NA se non si possono ottenere i valori
    } else {
      # --- 5. Calcolo dell'MSE per la run corrente ---
      bias_run <- stima_lambda - valore_obiettivo_lambda
      varianza_stimata_dal_modello_run <- se_lambda^2
      
      mse_lambda_runs[i] <- bias_run^2 + varianza_stimata_dal_modello_run
    }
  } else {
    mse_lambda_runs[i] <- NA # Segna come NA se il modello non ha fittato
  }
  
  # Stampa un aggiornamento periodico
  if (i %% 1 == 0 || i == n_simulazioni) {
    cat("Completata 'run' di fit:", i, "/", n_simulazioni, "\n")
  }
}

# --- 6. Calcolo dell'MSE medio ---
mse_medio_lambda_gmm <- mean(mse_lambda_runs, na.rm = TRUE)

cat("\n--- Risultati Finali ---\n")
print(paste("Valore obiettivo ipotetico per lambda:", valore_obiettivo_lambda))

# Stampa l'ultima stima valida di lambda e il suo SE (saranno rappresentative se i fit sono stabili)
if(exists("stima_lambda") && !is.null(stima_lambda)) {
  print(paste("Stima di Lambda (ultima run valida):", round(stima_lambda, 5)))
  print(paste("SE di Lambda (ultima run valida):", round(se_lambda, 5)))
}

print(paste("MSE medio calcolato per lambda (basato su SE del modello):", round(mse_medio_lambda_gmm, 7)))

# Nota: Se n_simulazioni > 1 e i dati di input non cambiano (come in questo script),
# tutti i valori in mse_lambda_runs dovrebbero essere identici, assumendo una convergenza deterministica.
# Per una vera simulazione Monte Carlo, dovresti generare nuovi dati 'y_gm' in ogni iterazione
# basati su parametri noti, o usare tecniche come il bootstrap sui dati originali.
# Questo script si concentra sull'illustrare il calcolo dell'MSE una volta che hai la stima e l'SE.
# print("Valori di MSE per ogni run (dovrebbero essere simili/uguali in questo esempio):")
# print(round(mse_lambda_runs, 7))

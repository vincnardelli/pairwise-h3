#-----------------------------------------------------------------------------#
# 0. Caricamento Librerie e Impostazioni Iniziali per DEBUG
#-----------------------------------------------------------------------------#
library(sf)
library(spdep)
library(spatialreg) # Non strettamente necessario per questo debug specifico di dati/plot
library(Matrix)
library(dplyr)
library(h3jsr)
library(ggplot2)
library(mvtnorm) # Per rmvnorm in generate_base_spatial_data

set.seed(12345) # Per riproducibilità del debug

#-----------------------------------------------------------------------------#
# A. DEFINIZIONI COMPLETE DELLE FUNZIONI HELPER
#-----------------------------------------------------------------------------#

#' Funzione per Generare u con Approssimazione Taylor/Neumann (efficiente)
generate_u_taylor_approx_efficient <- function(epsilon_vec, W_sparse_matrix, lambda, K_max, verbose = TRUE) {
  if (abs(lambda) < 1e-9 || K_max == 0) {
    if (verbose) cat("  Approssimazione Taylor (DGP): lambda=0 o K_max=0, u = epsilon.\n")
    return(epsilon_vec)
  }
  
  u_approx <- epsilon_vec 
  term_potenza_lambdaW_epsilon <- epsilon_vec 
  
  if (verbose) cat(paste0("  Approssimazione Taylor (DGP) K_max=", K_max, ": "))
  for (k_iter in 1:K_max) {
    if (verbose && (k_iter %% 10 == 0 || k_iter == K_max || K_max < 10 || k_iter == 1)) cat(paste0(k_iter, "..."))
    term_potenza_lambdaW_epsilon <- lambda * (W_sparse_matrix %*% term_potenza_lambdaW_epsilon)
    u_approx <- u_approx + term_potenza_lambdaW_epsilon
  }
  if (verbose) cat(" Fatto.\n")
  return(u_approx)
}

#' Genera dati spaziali di base 
generate_base_spatial_data <- function(n_points, data_type = "uniform", extent_polygon,
                                       x_mean = 0, x_sd = 1,
                                       cluster_lambda_centroids = 5,
                                       cluster_lambda_points_per_cluster = 100,
                                       cluster_sigma = 0.1,
                                       beta_true_vec, 
                                       lambda_true,   
                                       sigma_sq_eps_true, 
                                       k_neighbors_for_W,
                                       K_max_taylor_dgp, 
                                       verbose_dgp_approx = FALSE 
) {
  if (n_points == 0) {
    return(st_sf(geometry = st_sfc(crs = st_crs(extent_polygon)),
                 x_value = numeric(0), y_value = numeric(0), unique_id = integer(0)))
  }
  if (!inherits(extent_polygon, "sfc") && !inherits(extent_polygon, "sf")) {
    stop("extent_polygon deve essere un oggetto sfc o sf.")
  }
  if (inherits(extent_polygon, "sf")) {
    extent_polygon <- st_geometry(extent_polygon)
  }
  
  # --- 1. Generazione Punti ---
  points_sfc <- NULL
  if (data_type == "uniform") {
    points_sfc <- st_sample(extent_polygon, size = n_points, type = "random", exact = FALSE)
  } else if (data_type == "clustered") {
    num_centroids <- rpois(1, cluster_lambda_centroids)
    if (num_centroids == 0 && n_points > 0) num_centroids <- 1 # Assicura almeno un centroide se richiesti punti
    
    if (num_centroids == 0 ) { # Se n_points è 0 o num_centroids è 0 dopo rpois
      points_sfc <- st_sfc(crs = st_crs(extent_polygon)) 
    } else {
      centroid_points_sfc <- st_sample(extent_polygon, size = num_centroids, type = "random", exact = FALSE)
      all_cluster_points_list <- list()
      points_per_cluster_assigned <- rep(0, num_centroids)
      if (num_centroids > 0 && n_points > 0) { # Distribuisci punti solo se n_points > 0
        base_pts <- floor(n_points / num_centroids)
        remainder_pts <- n_points - (base_pts * num_centroids)
        points_per_cluster_assigned <- rep(base_pts, num_centroids)
        if (remainder_pts > 0) {
          safe_sample_indices <- if(num_centroids == 1) 1 else sample(num_centroids, remainder_pts)
          points_per_cluster_assigned[safe_sample_indices] <- points_per_cluster_assigned[safe_sample_indices] + 1
        }
      } # else points_per_cluster_assigned rimane zero se n_points è 0
      
      for (i in 1:num_centroids) {
        n_pts_this_cluster <- points_per_cluster_assigned[i]
        if (n_pts_this_cluster <= 0) next
        
        centroid_coord_matrix <- st_coordinates(centroid_points_sfc[i]) # Dovrebbe essere sempre matrice
        current_centroid_mean <- centroid_coord_matrix[1, , drop=FALSE] # Mantiene come matrice 1x2
        
        # DEBUG per l'errore 'names'
        cat(paste0("\nDEBUG ciclo cluster i=", i, ":\n"))
        cat(paste0("  n_pts_this_cluster = ", n_pts_this_cluster, "\n"))
        cat("  Classe di centroid_coord_matrix: ", class(centroid_coord_matrix),"; Dim: ", paste(dim(centroid_coord_matrix), collapse="x"), "\n")
        cat("  Valore di current_centroid_mean: ", paste(current_centroid_mean, collapse=", "), "\n")
        
        current_sigma_matrix <- diag(c(cluster_sigma^2, cluster_sigma^2))
        cat("  Dimensioni di current_sigma_matrix: ", paste(dim(current_sigma_matrix), collapse="x"), "\n")
        
        cluster_coords <- rmvnorm(n_pts_this_cluster,
                                  mean = current_centroid_mean[1,], # Prende il vettore riga
                                  sigma = current_sigma_matrix)
        
        cat("  Classe di cluster_coords: ", class(cluster_coords), "; Dim: ", paste(dim(cluster_coords), collapse="x"), "\n")
        if (is.matrix(cluster_coords)) {cat("  Prime righe di cluster_coords:\n"); print(head(cluster_coords))} else {print(cluster_coords)}
        
        points_df_this_cluster <- as.data.frame(cluster_coords)
        cat("  Classe di points_df_this_cluster: ", class(points_df_this_cluster), "; Dim: ", paste(dim(points_df_this_cluster), collapse="x"), "\n")
        cat("  Nomi colonne PRIMA: ", paste(names(points_df_this_cluster), collapse=", "), "\n")
        
        # ----> RIGA SOSPETTA PER L'ERRORE 'names' <----
        names(points_df_this_cluster) <- c("X", "Y") 
        cat("  Nomi colonne DOPO: ", paste(names(points_df_this_cluster), collapse=", "), "\n")
        # ----> FINE RIGA SOSPETTA <----
        
        cluster_sf_single <- st_as_sf(points_df_this_cluster, coords = c("X", "Y"), crs = st_crs(extent_polygon))
        all_cluster_points_list[[length(all_cluster_points_list) + 1]] <- cluster_sf_single
      } # Fine loop for
      
      if (length(all_cluster_points_list) > 0) {
        points_sfc_raw <- do.call(rbind, all_cluster_points_list) %>% st_geometry()
        valid_indices <- st_within(points_sfc_raw, extent_polygon, sparse = FALSE)[,1]
        points_sfc <- points_sfc_raw[valid_indices]
      } else if (n_points > 0) { 
        points_sfc <- st_sample(extent_polygon, size = n_points, type = "random", exact = FALSE)
      } else { 
        points_sfc <- st_sfc(crs = st_crs(extent_polygon)) 
      }
    }
  } else { stop("data_type non valido.") }
  
  current_n_final <- length(points_sfc)
  if (n_points > 0 && current_n_final == 0) { 
    points_sfc <- st_sample(extent_polygon, size = n_points, type = "random", exact = FALSE)
  } else if (current_n_final > n_points) { 
    points_sfc <- points_sfc[sample(current_n_final, n_points)]
  } else if (current_n_final < n_points && current_n_final >= 0) { 
    additional_needed <- n_points - current_n_final
    if (additional_needed > 0) { 
      fallback_pts <- st_sample(extent_polygon, size = additional_needed, type = "random", exact = FALSE)
      if(length(points_sfc) > 0 && length(fallback_pts) > 0) points_sfc <- c(points_sfc, fallback_pts)
      else if (length(fallback_pts) > 0) points_sfc <- fallback_pts
      # Se points_sfc era vuoto e fallback_pts è vuoto, rimane vuoto
    }
  }
  
  points_sf <- st_as_sf(points_sfc) # Questo può creare un sf object con 0 righe se points_sfc è vuoto
  actual_n_points <- nrow(points_sf)
  
  if (actual_n_points > 0) {
    geom_col_name <- attr(points_sf, "sf_column")
    if (!is.null(geom_col_name) && geom_col_name != "geometry" ) {
      names(points_sf)[names(points_sf) == geom_col_name] <- "geometry"
      st_geometry(points_sf) <- "geometry"
    } else if (is.null(geom_col_name) && "x" %in% names(points_sf) && inherits(points_sf$x, "sfc")) { 
      st_geometry(points_sf) <- "x"
      names(points_sf)[names(points_sf) == "x"] <- "geometry"
    } else if (is.null(geom_col_name) && !"geometry" %in% names(points_sf) && ncol(suppressWarnings(st_coordinates(points_sf))) > 0 && nrow(suppressWarnings(st_coordinates(points_sf))) == actual_n_points ) {
      warning("Colonna geometria non standard, la geometria potrebbe essere invalida.")
    }
    if (!"geometry" %in% names(points_sf) || (length(st_geometry(points_sf)) == 0 && actual_n_points > 0) ){
      warning(paste0("Geometria non valida o mancante per ", actual_n_points, " osservazioni in generate_base_spatial_data. Tentativo di return vuoto."))
      return(st_sf(geometry = st_sfc(crs = st_crs(extent_polygon)), x_value = numeric(0), y_value = numeric(0), unique_id = integer(0)))
    }
    
    points_sf$x_value <- rnorm(actual_n_points, mean = x_mean, sd = x_sd); points_sf$unique_id <- 1:actual_n_points
    lambda_eff = lambda_true; listw_global <- NULL
    if (actual_n_points > 1 && k_neighbors_for_W > 0) { 
      coords_matrix <- st_coordinates(points_sf)
      if(nrow(coords_matrix) == actual_n_points && actual_n_points > 0 && ncol(coords_matrix) == 2) {
        knn <- knearneigh(coords_matrix, k = min(k_neighbors_for_W, actual_n_points - 1)); nb <- knn2nb(knn)
        if(any(card(nb) == 0) && actual_n_points > 1) { warning("DGP: Alcuni punti non hanno vicini con k=", k_neighbors_for_W) }
        listw_global <- nb2listw(nb, style = "W", zero.policy = TRUE)
      } else if (actual_n_points > 0) {warning("Coordinate matrix row mismatch o no valid coords in DGP for W matrix creation.")}
    }
    X_matrix <- cbind(Intercept = 1, x_value = points_sf$x_value)
    if (ncol(X_matrix) != length(beta_true_vec)) {stop("Errore dimensioni beta_true_vec e X_matrix.")}
    mean_signal <- X_matrix %*% beta_true_vec
    epsilon <- rnorm(actual_n_points, mean = 0, sd = sqrt(sigma_sq_eps_true))
    if (!is.null(listw_global) && abs(lambda_eff) > 1e-9 && actual_n_points > 0 && K_max_taylor_dgp > 0) {
      W_mat_sparse <- as(listw_global, "CsparseMatrix")
      u <- generate_u_taylor_approx_efficient(epsilon, W_mat_sparse, lambda_eff, K_max_taylor_dgp, verbose = verbose_dgp_approx)
    } else { u <- epsilon }
    if (actual_n_points == 0) { points_sf$y_value <- numeric(0) } else { points_sf$y_value <- as.vector(mean_signal + u) }
  } else { 
    points_sf <- st_sf(geometry = st_sfc(crs = st_crs(extent_polygon)),
                       x_value = numeric(0), y_value = numeric(0), unique_id = integer(0))
  }
  return(points_sf)
}

#' Seleziona celle H3
select_h3_cells <- function(sf_data, h3_res, min_obs_in_h3,
                            k_ring_separation, target_selected_h3_count,
                            max_iterations_limit,
                            scenario_id_for_debug = "N/A",
                            debug_h3_selection = FALSE) {
  if (nrow(sf_data) == 0 || !"geometry" %in% names(sf_data) || length(st_geometry(sf_data)) == 0 ) {
    if(debug_h3_selection) cat(paste0("DEBUG H3 (", scenario_id_for_debug, "): Dati input vuoti o geometria mancante per select_h3_cells.\n"))
    return(character(0))
  }
  
  are_geoms_valid <- st_is_valid(sf_data)
  if(!all(are_geoms_valid)){
    if(debug_h3_selection) cat(paste0("WARN (", scenario_id_for_debug, "): ", sum(!are_geoms_valid), " geometrie non valide in sf_data. Tentativo di correzione con st_make_valid.\n"))
    sf_data <- st_make_valid(sf_data)
    if (any(!st_is_valid(sf_data))) { # Ricontrolla dopo la correzione
      cat(paste0("ERRORE (", scenario_id_for_debug, "): Correzione geometrie fallita. Impossibile procedere con select_h3_cells.\n"))
      return(character(0))
    }
  }
  
  sf_data_with_h3 <- sf_data %>%
    mutate(h3_index = as.character(h3jsr::point_to_cell(geometry, res = h3_res, simple=FALSE))) 
  
  h3_counts <- sf_data_with_h3 %>%
    st_drop_geometry() %>%
    group_by(h3_index) %>%
    summarise(n_obs = n(), .groups = 'drop') %>%
    filter(n_obs >= min_obs_in_h3)
  
  if (nrow(h3_counts) == 0) {
    if (debug_h3_selection) cat(paste0("DEBUG H3 (", scenario_id_for_debug, "): Nessuna cella H3 candidata iniziale con >= ", min_obs_in_h3, " oss.\n"))
    return(character(0))
  }
  candidate_h3s_master <- as.character(h3_counts$h3_index)
  selected_h3_list <- character(0)
  if (debug_h3_selection) {
    cat(paste0("\n--- DEBUG Inizio select_h3_cells (Scenario: ", scenario_id_for_debug, ") ---\n"))
    cat(paste0("      Numero candidate MASTRO iniziali (con >= ", min_obs_in_h3, " oss.): ", length(candidate_h3s_master), "\n"))
  }
  current_candidates_for_sampling <- candidate_h3s_master; iter_count_while <- 0
  while(length(selected_h3_list) < target_selected_h3_count &&
        length(current_candidates_for_sampling) > 0 &&
        iter_count_while < max_iterations_limit) {
    iter_count_while <- iter_count_while + 1
    if(debug_h3_selection && iter_count_while == 1) {cat(paste0("Prima del campionamento (iter 1): ", length(current_candidates_for_sampling), " candidati.\n"))}
    
    s_selected <- if (length(current_candidates_for_sampling) == 1) current_candidates_for_sampling[1] else sample(current_candidates_for_sampling, 1)
    selected_h3_list <- c(selected_h3_list, s_selected)
    
    cells_to_exclude_this_iteration <- as.character(s_selected)
    if (k_ring_separation > 0) {
      try_rings_call <- try(h3jsr::get_ring(h3_address = s_selected, ring_size = k_ring_separation, simple = TRUE), silent = TRUE)
      if (!inherits(try_rings_call, "try-error") && length(try_rings_call) > 0 && !is.null(unlist(try_rings_call))) { 
        cells_to_exclude_this_iteration <- unique(as.character(unlist(try_rings_call)))
      } else { 
        if (debug_h3_selection) {cat(paste0("DEBUG (", scenario_id_for_debug, "): Errore/output vuoto da get_ring per ", s_selected, ". Escludo solo s_selected.\n"))}
      }
    }
    current_candidates_for_sampling <- current_candidates_for_sampling[!(current_candidates_for_sampling %in% cells_to_exclude_this_iteration)]
    if (debug_h3_selection && (iter_count_while %% 10 == 0 || length(current_candidates_for_sampling) == 0 || length(selected_h3_list) == target_selected_h3_count) ) {cat(paste0("Iter ", iter_count_while, ", Selezionate ", length(selected_h3_list), ", Candidate rimaste ", length(current_candidates_for_sampling), "\n"))}
  }
  if (length(selected_h3_list) < target_selected_h3_count && iter_count_while < max_iterations_limit && length(current_candidates_for_sampling) == 0 ) {cat(paste0("WARN (", scenario_id_for_debug, "): Candidate esaurite. Selezionate: ", length(selected_h3_list), "/", target_selected_h3_count, "\n"))
  } else if (iter_count_while >= max_iterations_limit && length(selected_h3_list) < target_selected_h3_count) {cat(paste0("WARN (", scenario_id_for_debug, "): Limite iterazioni H3 raggiunto. Selezionate: ", length(selected_h3_list), "/", target_selected_h3_count, "\n"))}
  if (debug_h3_selection) {cat(paste0("\n--- DEBUG Fine select_h3_cells ---\nRestituite ", length(unique(selected_h3_list)), " celle H3 uniche.\n"))}
  return(unique(selected_h3_list)) 
}

#' Genera osservazioni accoppiate
generate_paired_observations <- function(all_sf_data, selected_h3_ids, h3_res) {
  if (length(selected_h3_ids) == 0) return(data.frame())
  if (!"y_value" %in% names(all_sf_data) || !"geometry" %in% names(all_sf_data) || length(st_geometry(all_sf_data)) == 0) {
    stop("Dati (y_value o geometry) mancanti o geometria vuota in all_sf_data per generate_paired_observations.")
  }
  
  # Assicura che st_is_valid sia chiamato prima di point_to_cell
  valid_geom_indices_gpo <- st_is_valid(all_sf_data$geometry, reason=FALSE)
  if(!all(valid_geom_indices_gpo)){
    # cat(paste0("WARN (generate_paired_observations): ", sum(!valid_geom_indices_gpo), " geometrie non valide. Filtraggio.\n"))
    all_sf_data <- all_sf_data[valid_geom_indices_gpo, ]
    if(nrow(all_sf_data) == 0) return(data.frame())
  }
  
  temp_all_sf_data_with_h3 <- all_sf_data %>%
    mutate(h3_index_for_pairing = as.character(h3jsr::point_to_cell(geometry, res = h3_res, simple=FALSE)))
  
  paired_observations_list <- list()
  unique_selected_h3_ids <- unique(selected_h3_ids) # Assicura di iterare solo su ID unici
  
  for (h3_id_current in unique_selected_h3_ids) {
    points_in_h3 <- temp_all_sf_data_with_h3 %>% 
      filter(h3_index_for_pairing == as.character(h3_id_current))
    
    if (nrow(points_in_h3) >= 2) {
      pair_indices <- sample(nrow(points_in_h3), 2, replace = FALSE)
      obs_i <- points_in_h3[pair_indices[1], ]; obs_l <- points_in_h3[pair_indices[2], ]
      
      paired_observations_list[[length(paired_observations_list) + 1]] <- data.frame(
        h3_id_selected = h3_id_current, 
        unique_id_i = obs_i$unique_id, x_i = obs_i$x_value, y_i = obs_i$y_value,
        unique_id_l = obs_l$unique_id, x_l = obs_l$x_value, y_l = obs_l$y_value)
    }
  }
  if (length(paired_observations_list) > 0) return(do.call(rbind, paired_observations_list))
  else return(data.frame())
}


#-----------------------------------------------------------------------------#
# B. PARAMETRI SPECIFICI PER QUESTO SCRIPT DI DEBUG
#-----------------------------------------------------------------------------#
N_debug <- 5000         # Ridotto per debug più veloce
H3_RES_debug <- 6        # RISOLUZIONE PIU' BASSA (esagoni più grandi)
LAMBDA_DGP_debug <- -0.3 
K_MAX_TAYLOR_debug <- 30 
TARGET_Q_PAIRS_debug <- 50 
MIN_OBS_H3_debug <- 2
K_RING_SEP_debug <- 1    
K_NEIGH_W_dgp_debug <- 6  
BETA_VEC_debug <- c(1.0, 1.5) # Intercetta, Beta1
SIGMA_SQ_EPS_DGP_debug <- 1.0

cat("--- INIZIO SCRIPT DI DEBUG ---\n")
cat(paste("Parametri Debug: N =", N_debug, ", H3_Res =", H3_RES_debug, 
          ", Lambda_DGP =", LAMBDA_DGP_debug, ", K_Taylor =", K_MAX_TAYLOR_debug,
          ", TargetQ =", TARGET_Q_PAIRS_debug, ", K_Ring_Sep =", K_RING_SEP_debug, "\n\n"))

#-----------------------------------------------------------------------------#
# C. ESECUZIONE PASSO-PASSO DELLE FASI CHIAVE
#-----------------------------------------------------------------------------#

# 1. Genera dati base
cat("--- FASE 1: Generazione Dati Base ---\n")
debug_extent_polygon <- st_polygon(list(matrix(c(0,0, 1,0, 1,1, 0,1, 0,0), ncol=2, byrow=TRUE))) %>% st_sfc(crs = 4326)

base_data_debug <- generate_base_spatial_data(
  n_points = N_debug,
  data_type = "uniform", 
  extent_polygon = debug_extent_polygon,
  beta_true_vec = BETA_VEC_debug,
  lambda_true = LAMBDA_DGP_debug,
  sigma_sq_eps_true = SIGMA_SQ_EPS_DGP_debug,
  k_neighbors_for_W = K_NEIGH_W_dgp_debug,
  K_max_taylor_dgp = K_MAX_TAYLOR_debug,
  verbose_dgp_approx = TRUE 
)
cat(paste("Dati base generati:", nrow(base_data_debug), "osservazioni.\n"))
if(nrow(base_data_debug) == 0 && N_debug > 0) stop("Generazione dati base fallita, 0 osservazioni prodotte.")
print(head(select(base_data_debug, unique_id, x_value, y_value, geometry)))
cat("\n")

# 2. Seleziona celle H3
cat("--- FASE 2: Selezione Celle H3 ---\n")
selected_h3_ids_debug <- select_h3_cells(
  sf_data = base_data_debug,
  h3_res = H3_RES_debug,
  min_obs_in_h3 = MIN_OBS_H3_debug,
  k_ring_separation = K_RING_SEP_debug,
  target_selected_h3_count = TARGET_Q_PAIRS_debug,
  max_iterations_limit = N_debug, 
  scenario_id_for_debug = "DEBUG_SCRIPT",
  debug_h3_selection = TRUE 
)
cat(paste("Numero di celle H3 selezionate (output da select_h3_cells):", length(selected_h3_ids_debug), "\n"))
cat(paste("Numero di ID H3 UNICI selezionati:", length(unique(selected_h3_ids_debug)), "\n"))
if (length(selected_h3_ids_debug) != length(unique(selected_h3_ids_debug))) {
  cat("ATTENZIONE: Ci sono ID H3 duplicati in selected_h3_ids_debug!\n")
  print(paste("ID duplicati:", names(table(selected_h3_ids_debug)[table(selected_h3_ids_debug)>1])))
} else {
  cat("Conferma: Tutti gli ID H3 selezionati da select_h3_cells sono unici (grazie a return unique()).\n")
}
cat("Primi 10 H3 selezionati:", head(selected_h3_ids_debug, 10), "\n\n")


if (length(selected_h3_ids_debug) == 0) {
  message("Nessuna cella H3 selezionata. Impossibile procedere con il debug delle coppie e del plot.")
} else {
  # 3. Genera osservazioni accoppiate
  cat("--- FASE 3: Generazione Osservazioni Accoppiate ---\n")
  paired_obs_df_debug <- generate_paired_observations(
    all_sf_data = base_data_debug,
    selected_h3_ids = selected_h3_ids_debug, # select_h3_cells ora dovrebbe restituire unici
    h3_res = H3_RES_debug
  )
  cat(paste("Numero di coppie generate (righe in paired_obs_df_debug):", nrow(paired_obs_df_debug), "\n"))
  if (nrow(paired_obs_df_debug) > 0) {
    cat(paste("Numero di h3_id_selected UNICI nel dataframe delle coppie:", 
              length(unique(paired_obs_df_debug$h3_id_selected)), "\n"))
    
    # Dato che generate_paired_observations ora itera su unique(selected_h3_ids) internamente,
    # questo controllo verifica se il numero di coppie corrisponde al numero di H3 uniche idonee.
    if (length(unique(selected_h3_ids_debug)) >= length(unique(paired_obs_df_debug$h3_id_selected))) {
      cat("Conferma: Il numero di H3 uniche nel df delle coppie è <= al numero di H3 uniche selezionate (alcune H3 potrebbero non aver avuto >=2 punti).\n")
    } else {
      cat("ATTENZIONE: Più H3 uniche nel df delle coppie che quelle selezionate! Improbabile.\n")
    }
    count_coppie_per_h3 <- paired_obs_df_debug %>% group_by(h3_id_selected) %>% summarise(n_coppie = n())
    if(any(count_coppie_per_h3$n_coppie > 1)) {
      cat("ATTENZIONE GRAVE: Alcune celle H3 hanno generato più di una coppia! Questo non dovrebbe succedere se selected_h3_ids è unico E generate_paired_observations itera su ID unici.\n")
      print(filter(count_coppie_per_h3, n_coppie > 1))
    } else {
      cat("Conferma: Ogni cella H3 unica ha generato al massimo una coppia.\n")
    }
  } else {
    cat("Nessuna coppia generata.\n")
  }
  cat("\n")
  
  # 4. Identifica i punti "rossi" per il plot
  cat("--- FASE 4: Identificazione Punti per il Plot ---\n")
  paired_points_ids_debug <- character(0)
  if (nrow(paired_obs_df_debug) > 0) {
    paired_points_ids_debug <- unique(c(paired_obs_df_debug$unique_id_i, paired_obs_df_debug$unique_id_l))
  }
  cat(paste("Numero di unique_id di punti usati nelle coppie (punti candidati ad essere rossi):", length(paired_points_ids_debug), "\n"))
  
  paired_points_sf_debug <- base_data_debug %>% 
    filter(unique_id %in% paired_points_ids_debug)
  cat(paste("Numero di righe (punti) in paired_points_sf_debug:", nrow(paired_points_sf_debug), "\n"))
  cat("\n")
  
  #-----------------------------------------------------------------------------#
  # D. ANALISI SPECIFICA DI UNA CELLA H3 "PROBLEMATICA"
  #-----------------------------------------------------------------------------#
  cat("--- FASE 5: Analisi Dettagliata Celle e Plot ---\n")
  
  ID_CELLA_PER_PLOT_DETTAGLIATO <- NULL
  # Se vuoi forzare una cella specifica dal tuo grafico originale, mettila qui:
  # ID_CELLA_PER_PLOT_DETTAGLIATO <- "86754e807ffffff" # Esempio dal tuo titolo, risoluzione 6
  
  if (is.null(ID_CELLA_PER_PLOT_DETTAGLIATO) && length(selected_h3_ids_debug) > 0) {
    # Prova a prendere la prima cella H3 selezionata per l'analisi dettagliata
    # o una che ha molti punti, se possibile identificarla.
    # Qui prendiamo la prima per semplicità.
    ID_CELLA_PER_PLOT_DETTAGLIATO <- selected_h3_ids_debug[1] 
    cat(paste("Cella H3 scelta automaticamente per il plot dettagliato:", ID_CELLA_PER_PLOT_DETTAGLIATO, "\n"))
  } else if (length(selected_h3_ids_debug) == 0) {
    cat("Nessuna cella H3 selezionata, impossibile fare plot dettagliato.\n")
  }
  
  if (!is.null(ID_CELLA_PER_PLOT_DETTAGLIATO)) {
    ID_CELLA_PULITO <- trimws(ID_CELLA_PER_PLOT_DETTAGLIATO)
    is_valid_h3_check <- all(h3jsr::is_valid(ID_CELLA_PULITO)) # is_valid può restituire un vettore
    cat(paste0("L'ID H3 '", ID_CELLA_PULITO, "' è valido? : ", is_valid_h3_check, "\n"))
    
    if(!is_valid_h3_check){
      cat("ATTENZIONE: L'ID cella sospetta non è valido. Scegline uno valido da 'selected_h3_ids_debug'.\n")
    } else {
      base_data_debug_with_h3_assigned <- base_data_debug %>%
        mutate(h3_viz_debug = as.character(h3jsr::point_to_cell(geometry, res = H3_RES_debug, simple=FALSE)))
      
      punti_base_in_cella_target <- base_data_debug_with_h3_assigned %>%
        filter(h3_viz_debug == ID_CELLA_PULITO)
      cat(paste("Numero totale di punti base in cella", ID_CELLA_PULITO, ":", nrow(punti_base_in_cella_target), "\n"))
      cat("  Unique IDs dei punti base in cella target (max 10 mostrati):", paste(head(punti_base_in_cella_target$unique_id,10), collapse=", "), if(nrow(punti_base_in_cella_target)>10) "...", "\n")
      
      coppia_campionata_da_questa_cella <- paired_obs_df_debug %>%
        filter(h3_id_selected == ID_CELLA_PULITO)
      cat(paste("Numero di coppie campionate da cella", ID_CELLA_PULITO, ":", nrow(coppia_campionata_da_questa_cella), "(dovrebbe essere 1 se selezionata e idonea)\n"))
      
      punti_rossi_DA_QUESTA_COPPIA_ID <- c()
      if (nrow(coppia_campionata_da_questa_cella) == 1) {
        punti_rossi_DA_QUESTA_COPPIA_ID <- c(coppia_campionata_da_questa_cella$unique_id_i, coppia_campionata_da_questa_cella$unique_id_l)
        cat("  Unique IDs dei 2 punti campionati per la coppia da questa cella:", 
            paste(punti_rossi_DA_QUESTA_COPPIA_ID, collapse=", "), "\n")
      }
      
      cella_target_poly <- NULL
      cella_target_poly <- cell_to_polygon(ID_CELLA_PULITO, simple = FALSE)
      
      punti_rossi_GEOMETRICAMENTE_in_cella_target <- st_sf(data.frame(), geometry = st_sfc(crs=st_crs(base_data_debug))) 
      if (nrow(paired_points_sf_debug) > 0 && !is.null(cella_target_poly) && nrow(cella_target_poly) > 0 && !st_is_empty(st_geometry(cella_target_poly))) {
        punti_rossi_GEOMETRICAMENTE_in_cella_target <- st_filter(paired_points_sf_debug, cella_target_poly, .predicate = st_intersects)
      }
      cat(paste("Numero di punti rossi (da tutti i punti accoppiati) che INTERSECANO geometricamente il poligono di", 
                ID_CELLA_PULITO, ":", nrow(punti_rossi_GEOMETRICAMENTE_in_cella_target), "\n"))
      if (nrow(punti_rossi_GEOMETRICAMENTE_in_cella_target) > 0) {
        cat("  Unique IDs dei punti rossi INTERSECANTI:", 
            paste(punti_rossi_GEOMETRICAMENTE_in_cella_target$unique_id, collapse=", "), "\n")
      }
      
      cat("Creazione plot dettagliato...\n")
      vicini_sospetta_poly <- NULL
      tryCatch({
        vicini_ids_list <- h3jsr::get_ring(ID_CELLA_PULITO, ring_size = 1) 
        if(length(vicini_ids_list[[1]]) > 0) vicini_sospetta_poly <- cell_to_polygon(unlist(vicini_ids_list), simple = FALSE)
      }, error = function(e) { cat("Nota: impossibile ottenere vicini per il plot.\n")})
      
      bbox_cella <- st_bbox(cella_target_poly)
      # Calcola il centro e l'estensione della cella per lo zoom
      # Adatta lo span in base alla risoluzione H3 per uno zoom appropriato
      span_factor_res <- case_when(H3_RES_debug <= 4 ~ 0.2, H3_RES_debug == 5 ~ 0.05, H3_RES_debug == 6 ~ 0.01, TRUE ~ 0.005)
      span_x_calc <- (bbox_cella["xmax"] - bbox_cella["xmin"])
      span_y_calc <- (bbox_cella["ymax"] - bbox_cella["ymin"])
      
      # Assicura che lo span non sia zero se la cella è degenere (improbabile per H3)
      if(span_x_calc < 1e-6) span_x_calc <- span_factor_res * 5
      if(span_y_calc < 1e-6) span_y_calc <- span_factor_res * 5
      
      plot_xlim <- c(bbox_cella["xmin"] - span_x_calc*1.5, bbox_cella["xmax"] + span_x_calc*1.5)
      plot_ylim <- c(bbox_cella["ymin"] - span_y_calc*1.5, bbox_cella["ymax"] + span_y_calc*1.5)
      
      plot_bbox_sf <- st_as_sfc(st_bbox(c(xmin=plot_xlim[1], ymin=plot_ylim[1], 
                                          xmax=plot_xlim[2], ymax=plot_ylim[2]), 
                                        crs=st_crs(base_data_debug)))
      
      base_data_plot_subset <- st_filter(base_data_debug_with_h3_assigned, plot_bbox_sf, .predicate = st_intersects)
      paired_points_plot_subset <- st_filter(paired_points_sf_debug, plot_bbox_sf, .predicate = st_intersects)
      
      p_debug <- ggplot()
      if(!is.null(vicini_sospetta_poly) && nrow(vicini_sospetta_poly) > 0) {
        p_debug <- p_debug + geom_sf(data = vicini_sospetta_poly, fill = NA, color = "grey70", linetype = "dotted", linewidth=0.3) # Era size
      }
      p_debug <- p_debug +
        geom_sf(data = base_data_plot_subset, aes(shape = "Tutti i Punti Base"), color="grey60", size = 1.5, alpha = 0.4) +
        geom_sf(data = cella_target_poly, aes(fill = "Cella Target Analizzata"), alpha = 0.15, color = "blue", linewidth=1.0) # Era size
      
      if (nrow(paired_points_plot_subset) > 0) {
        p_debug <- p_debug + 
          geom_sf(data = paired_points_plot_subset, aes(shape = "Punti nelle Coppie (Rossi)"), color = "red", size = 3, stroke=1)
      }
      
      if (length(punti_rossi_DA_QUESTA_COPPIA_ID) > 0) { 
        punti_specifici_coppia <- base_data_debug %>% filter(unique_id %in% punti_rossi_DA_QUESTA_COPPIA_ID)
        if (nrow(punti_specifici_coppia) > 0) {
          p_debug <- p_debug + 
            geom_sf(data = punti_specifici_coppia, aes(shape="Coppia da Cella Target"), color="magenta", size=4, stroke=1.5)
        }
      }
      
      p_debug <- p_debug + 
        scale_shape_manual(name="Legenda Punti", 
                           values=c("Tutti i Punti Base"=1, "Punti nelle Coppie (Rossi)"=16, "Coppia da Cella Target"=8), # 8 è una stella
                           labels=c("Tutti i Punti Base", "Punti in una Coppia (da qualsiasi cella)", "Coppia Specifica da Cella Target")) +
        scale_fill_manual(name="Legenda Celle", values=c("Cella Target Analizzata"="orange")) +
        coord_sf(xlim = plot_xlim, ylim = plot_ylim, expand = TRUE) + 
        labs(title = paste("Debug Dettagliato Cella H3:", ID_CELLA_PULITO, "(Res:", H3_RES_debug,")"),
             subtitle = paste(nrow(punti_base_in_cella_target), "punti base in cella.",
                              nrow(coppia_campionata_da_questa_cella),"coppia camp. da questa cella ->", length(punti_rossi_DA_QUESTA_COPPIA_ID), "punti magenta.",
                              nrow(punti_rossi_GEOMETRICAMENTE_in_cella_target), "punti rossi totali intersecano geometricamente questa cella."),
             x = "Longitudine", y = "Latitudine") +
        theme_minimal(base_size = 10) +
        theme(legend.position = "bottom", legend.box="vertical", plot.title = element_text(size=10), plot.subtitle = element_text(size=8))
      
      print(p_debug)
      
      debug_plot_filename <- paste0("debug_plot_H3_", gsub("[:/]", "_", ID_CELLA_PULITO), 
                                    "_N", N_debug, "_Res", H3_RES_debug, ".png")
      ggsave(file.path(tempdir(), debug_plot_filename), plot = p_debug, width = 10, height = 10, dpi=150)
      cat(paste0("Plot di debug salvato in: ", file.path(tempdir(), debug_plot_filename), "\n"))
    }
  } else {
    cat("ID_CELLA_PER_PLOT_DETTAGLIATO non impostato o non valido. Impossibile creare plot dettagliato.\n")
  }
} # Fine del blocco 'else' per selected_h3_ids_debug > 0

cat("\n--- FINE SCRIPT DI DEBUG ---\n")
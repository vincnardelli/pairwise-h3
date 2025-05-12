library(digest)
library(readr)
library(sf)
library(dplyr)
library(h3jsr)
library(mvtnorm)
library(ggplot2)
library(tidyr)
library(progress)
library(Matrix)
library(spatialreg) 
library(spdep)      

#-----------------------------------------------------------------------------#
# 1. IMPOSTAZIONI GLOBALI E PARAMETRI DI SIMULAZIONE ----
# (Questi verranno definiti dopo le funzioni helper)
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# 2. FUNZIONI HELPER ----
#-----------------------------------------------------------------------------#

#' Genera dati spaziali di base (posizioni e covariata x)
generate_base_spatial_data <- function(n_points, data_type = "uniform", extent_polygon,
                                       x_mean = 0, x_sd = 1,
                                       cluster_lambda_centroids = 5,
                                       cluster_lambda_points_per_cluster = 100,
                                       cluster_sigma = 0.1) {
  if (n_points == 0) { # Gestione esplicita all'inizio
    return(st_sf(geometry = st_sfc(crs = st_crs(extent_polygon)),
                 x_value = numeric(0),
                 unique_id = integer(0)))
  }
  if (!inherits(extent_polygon, "sfc") && !inherits(extent_polygon, "sf")) {
    stop("extent_polygon deve essere un oggetto sfc o sf.")
  }
  if (inherits(extent_polygon, "sf")) {
    extent_polygon <- st_geometry(extent_polygon)
  }
  
  if (data_type == "uniform") {
    points_sfc <- st_sample(extent_polygon, size = n_points, type = "random")
  } else if (data_type == "clustered") {
    num_centroids <- rpois(1, cluster_lambda_centroids)
    if (num_centroids == 0 && n_points > 0) num_centroids <- 1
    
    if (num_centroids == 0 ) { # Se n_points > 0, num_centroids è già 1. Questo copre n_points = 0 se non gestito sopra.
      points_sfc <- st_sfc(crs = st_crs(extent_polygon)) # Vuoto ma con CRS
    } else {
      centroid_points_sfc <- st_sample(extent_polygon, size = num_centroids, type = "random")
      all_cluster_points_list <- list()
      points_per_cluster_assigned <- rep(0, num_centroids)
      
      if (num_centroids > 0) { # Dovrebbe essere sempre vero qui
        base_pts <- floor(n_points / num_centroids)
        remainder_pts <- n_points - (base_pts * num_centroids)
        points_per_cluster_assigned <- rep(base_pts, num_centroids)
        if (remainder_pts > 0) {
          safe_sample_indices <- if(num_centroids == 1) 1 else sample(num_centroids, remainder_pts)
          points_per_cluster_assigned[safe_sample_indices] <- points_per_cluster_assigned[safe_sample_indices] + 1
        }
      }
      
      for (i in 1:num_centroids) {
        n_pts_this_cluster <- points_per_cluster_assigned[i]
        if (n_pts_this_cluster <= 0) next
        
        centroid_coord <- st_coordinates(centroid_points_sfc[i])
        cluster_coords <- rmvnorm(n_pts_this_cluster,
                                  mean = centroid_coord,
                                  sigma = diag(c(cluster_sigma^2, cluster_sigma^2)))
        
        points_df_this_cluster <- as.data.frame(cluster_coords)
        names(points_df_this_cluster) <- c("X", "Y")
        
        cluster_sf_single <- st_as_sf(points_df_this_cluster, coords = c("X", "Y"), crs = st_crs(extent_polygon))
        all_cluster_points_list[[length(all_cluster_points_list) + 1]] <- cluster_sf_single
      }
      
      if (length(all_cluster_points_list) > 0) {
        points_sfc_raw <- do.call(rbind, all_cluster_points_list) %>% st_geometry()
        points_sfc <- points_sfc_raw[st_intersects(points_sfc_raw, extent_polygon, sparse = FALSE)[,1]]
        
        current_n <- length(points_sfc)
        if (current_n > n_points) {
          points_sfc <- points_sfc[sample(current_n, n_points)]
        } else if (current_n < n_points && current_n >= 0) { 
          additional_needed <- n_points - current_n
          if (additional_needed > 0) {
            fallback_pts <- st_sample(extent_polygon, size = additional_needed, type = "random")
            points_sfc <- c(points_sfc, fallback_pts)
          }
        } 
      } else if (n_points > 0) { 
        points_sfc <- st_sample(extent_polygon, size = n_points, type = "random")
      } else { 
        points_sfc <- st_sfc(crs = st_crs(extent_polygon)) 
      }
    }
  } else {
    stop("data_type non valido. Scegliere 'uniform' o 'clustered'.")
  }
  
  current_n_final <- length(points_sfc)
  if (n_points > 0 && current_n_final == 0) { # Se n_points > 0 ma non sono stati generati punti
    points_sfc <- st_sample(extent_polygon, size = n_points, type = "random")
  } else if (current_n_final > n_points) {
    points_sfc <- points_sfc[sample(current_n_final, n_points)]
  } else if (current_n_final < n_points) { 
    additional_needed <- n_points - current_n_final
    if (additional_needed > 0) {
      fallback_pts <- st_sample(extent_polygon, size = additional_needed, type = "random")
      points_sfc <- c(points_sfc, fallback_pts)
    }
  }
  
  points_sf <- st_as_sf(points_sfc) # Può essere vuoto se n_points = 0
  if (nrow(points_sf) > 0) {
    geom_col_name <- attr(points_sf, "sf_column")
    if (geom_col_name != "geometry") {
      names(points_sf)[names(points_sf) == geom_col_name] <- "geometry"
      st_geometry(points_sf) <- "geometry"
    }
    points_sf$x_value <- rnorm(nrow(points_sf), mean = x_mean, sd = x_sd)
    points_sf$unique_id <- 1:nrow(points_sf)
  } else { # Assicura la struttura corretta se vuoto
    points_sf <- st_sf(geometry = st_sfc(crs = st_crs(extent_polygon)),
                       x_value = numeric(0),
                       unique_id = integer(0))
  }
  return(points_sf)
}


#' Seleziona celle H3
select_h3_cells <- function(sf_data, h3_res, min_obs_in_h3,
                            k_ring_separation, target_selected_h3_count,
                            max_iterations_limit,
                            scenario_id_for_debug = "N/A",
                            debug_h3_selection = FALSE) { 
  
  if (nrow(sf_data) == 0) return(character(0))
  
  sf_data_with_h3 <- sf_data %>%
    mutate(h3_index = as.character(h3jsr::point_to_cell(geometry, res = h3_res))) 
  
  h3_counts <- sf_data_with_h3 %>%
    st_drop_geometry() %>%
    group_by(h3_index) %>%
    summarise(n_obs = n(), .groups = 'drop') %>%
    filter(n_obs >= min_obs_in_h3)
  
  if (nrow(h3_counts) == 0) {
    if (debug_h3_selection) cat(paste0("DEBUG H3 (", scenario_id_for_debug, "): Nessuna cella H3 candidata iniziale con >= ", min_obs_in_h3, " oss.\n"))
    return(character(0))
  }
  
  # candidate_h3s_master contiene tutte le celle che soddisfano min_obs_in_h3
  candidate_h3s_master <- as.character(h3_counts$h3_index) 
  selected_h3_list <- character(0)
  
  if (debug_h3_selection) {
    cat(paste0("\n--- DEBUG Inizio select_h3_cells (Scenario: ", scenario_id_for_debug, ") ---\n"))
    cat(paste0("    Numero candidate MASTRO iniziali (con >= ", min_obs_in_h3, " oss.): ", length(candidate_h3s_master), "\n"))
  }
  
  # current_candidates_for_sampling è il pool da cui si campiona, viene ridotto ad ogni iterazione
  current_candidates_for_sampling <- candidate_h3s_master
  
  iter_count_while <- 0
  
  while(length(selected_h3_list) < target_selected_h3_count &&
        length(current_candidates_for_sampling) > 0 && 
        iter_count_while < max_iterations_limit) {
    
    iter_count_while <- iter_count_while + 1
    
    # Salva lo snapshot di current_candidates_for_sampling per il debug *prima* del campionamento
    if(debug_h3_selection) current_candidates_snapshot_before_sample <- current_candidates_for_sampling
    
    s_selected <- sample(current_candidates_for_sampling, 1) 
    selected_h3_list <- c(selected_h3_list, s_selected)
    
    cells_to_exclude_this_iteration <- as.character(s_selected) 
    
    if (k_ring_separation >= 0) { # k_ring_separation=0 significa escludere solo s_selected
      if (k_ring_separation > 0) {
        try_rings_call <- try(
          h3jsr::get_ring(h3_address = s_selected, ring_size = k_ring_separation, simple = TRUE),
          silent = TRUE # Mantenere silent=TRUE per evitare di fermare lo script, gestire l'errore dopo
        )
        if (!inherits(try_rings_call, "try-error") && length(try_rings_call) > 0) {
          cells_to_exclude_this_iteration <- unique(as.character(unlist(try_rings_call))) 
        } else {
          if (debug_h3_selection) {
            cat(paste0("        DEBUG (", scenario_id_for_debug, "): Errore o output vuoto da get_ring per s_selected=", s_selected," (k=",k_ring_separation,"). Escludo solo s_selected.\n"))
            if(inherits(try_rings_call, "try-error")) print(attributes(try_rings_call)$condition$message)
          }
        }
      }
    }
    
    if (debug_h3_selection) {
      cat(paste0("\n    --- DEBUG Iterazione While Loop: ", iter_count_while, " (Target: ",target_selected_h3_count ,", Selezionate: ", length(selected_h3_list),") ---\n"))
      cat(paste0("        Candidate disponibili per campionamento (current_candidates_snapshot_before_sample): ", length(current_candidates_snapshot_before_sample), "\n"))
      cat(paste0("        Cella Selezionata (s_selected): ", s_selected, "\n"))
      cat(paste0("        Celle da escludere (s_selected + k-ring k=",k_ring_separation,"): ", length(cells_to_exclude_this_iteration), " celle. Prime 10: ", paste(head(cells_to_exclude_this_iteration,10), collapse=", "), "\n"))
      
      # Celle che sono in cells_to_exclude E erano nel pool di campionamento di questa iterazione
      common_to_remove_now <- intersect(cells_to_exclude_this_iteration, current_candidates_snapshot_before_sample)
      cat(paste0("        Celle comuni tra 'cells_to_exclude' e 'current_candidates_snapshot_before_sample' (effettivamente rimosse da questo pool): ", length(common_to_remove_now), "\n"))
      if(length(common_to_remove_now) > 0 && debug_h3_selection) print(head(common_to_remove_now))
    }
    
    original_length_current_candidates <- length(current_candidates_for_sampling)
    # Rimuovi le celle da escludere dal set da cui si campionerà nelle prossime iterazioni
    current_candidates_for_sampling <- current_candidates_for_sampling[!(current_candidates_for_sampling %in% cells_to_exclude_this_iteration)]
    new_length_current_candidates <- length(current_candidates_for_sampling)
    
    if (debug_h3_selection) {
      cat(paste0("        DEBUG current_candidates_for_sampling dopo esclusione: ", new_length_current_candidates, " (riduzione di ", original_length_current_candidates - new_length_current_candidates ,")\n"))
      if ( (original_length_current_candidates - new_length_current_candidates) == 0 && 
           length(intersect(cells_to_exclude_this_iteration, current_candidates_snapshot_before_sample)) > 0 &&
           original_length_current_candidates > 0) { # Se non è diminuito ma doveva
        cat("        ---> DEBUG ALERT: Nessuna cella rimossa da current_candidates_for_sampling MA C'ERANO CELLE COMUNI DA RIMUOVERE!\n")
        cat("        ---> s_selected: ", s_selected, "\n")
      }
    }
  } # Fine while
  
  # ... (resto della funzione come prima, con i messaggi WARN)
  if (length(selected_h3_list) < target_selected_h3_count && length(current_candidates_for_sampling) == 0 && iter_count_while < max_iterations_limit) {
    cat(paste0("WARN (", scenario_id_for_debug, "): Candidate esaurite prima di raggiungere target_selected_h3_count. Selezionate: ", length(selected_h3_list), "/", target_selected_h3_count, "\n"))
  } else if (iter_count_while >= max_iterations_limit && length(selected_h3_list) < target_selected_h3_count) {
    cat(paste0("WARN (", scenario_id_for_debug, "): Raggiunto limite max_iterations_limit (", max_iterations_limit, ") prima di raggiungere target_selected_h3_count. Selezionate: ", length(selected_h3_list), "/", target_selected_h3_count, "\n"))
  }
  
  if (debug_h3_selection) {
    cat(paste0("\n--- DEBUG Fine select_h3_cells (Scenario: ", scenario_id_for_debug, ") ---\nSelezionate ", length(selected_h3_list), " celle H3 finali.\n"))
    if (length(selected_h3_list) > 0) {
      if (length(selected_h3_list) <= 50) print(selected_h3_list) else print(paste0("Prime 50 selezionate: ", paste(head(selected_h3_list, 50), collapse=", ")))
    }
    cat("--------------------------------\n\n")
  }
  return(selected_h3_list)
}

#' Genera osservazioni accoppiate e calcola la variabile dipendente y_paired
generate_paired_observations <- function(all_sf_data, selected_h3_ids, h3_res,
                                         beta_true, sigma_sq_true, psi_true) {
  if (length(selected_h3_ids) == 0) {
    return(data.frame())
  }
  
  # Assicura che h3_index sia presente e sia carattere
  if (!"h3_index" %in% names(all_sf_data)) {
    all_sf_data_with_h3 <- all_sf_data %>%
      mutate(h3_index = as.character(h3jsr::point_to_cell(geometry, res = h3_res)))
  } else {
    all_sf_data_with_h3 <- all_sf_data %>% mutate(h3_index = as.character(h3_index)) 
  }
  
  paired_observations_list <- list()
  
  for (h3_id_current in selected_h3_ids) {
    # Filtra i punti per la cella H3 corrente (assicurati che il confronto sia tra caratteri)
    points_in_h3 <- all_sf_data_with_h3 %>% filter(h3_index == as.character(h3_id_current))
    
    if (nrow(points_in_h3) >= 2) {
      pair_indices <- sample(nrow(points_in_h3), 2, replace = FALSE)
      obs_i <- points_in_h3[pair_indices[1], ]
      obs_l <- points_in_h3[pair_indices[2], ]
      
      cov_matrix_errors <- matrix(c(sigma_sq_true, psi_true * sigma_sq_true,
                                    psi_true * sigma_sq_true, sigma_sq_true), nrow = 2)
      
      errors_il <- tryCatch({
        rmvnorm(n = 1, mean = c(0,0), sigma = cov_matrix_errors)
      }, error = function(e) {
        warning(paste0("Errore in rmvnorm per H3 ", h3_id_current,
                       " con psi_true = ", psi_true, ", sigma_sq_true = ", sigma_sq_true,
                       ". Matrice di covarianza potrebbe non essere definita positiva. Uso errori non correlati. Messaggio: ", e$message))
        return(rnorm(2, mean = 0, sd = sqrt(sigma_sq_true))) # Fallback
      })
      
      e_i <- errors_il[1]
      e_l <- errors_il[2]
      
      # Calcola y_i e y_l: y = beta*x + errore_correlato_e + rumore_iid_N(0,1)
      # Il termine rnorm(1,0,1) è un rumore addizionale i.i.d. con varianza 1.
      # Questo rumore è SEPARATO da sigma_sq_true (che è la varianza dei termini e_i ed e_l).
      y_i <- beta_true * obs_i$x_value + e_i + rnorm(n = 1, mean = 0, sd = 1)
      y_l <- beta_true * obs_l$x_value + e_l + rnorm(n = 1, mean = 0, sd = 1)
      
      paired_observations_list[[length(paired_observations_list) + 1]] <- data.frame(
        h3_id_selected = h3_id_current,
        unique_id_i = obs_i$unique_id,
        x_i = obs_i$x_value,
        y_i = y_i,
        unique_id_l = obs_l$unique_id,
        x_l = obs_l$x_value,
        y_l = y_l
      )
    }
  }
  
  if (length(paired_observations_list) > 0) {
    return(do.call(rbind, paired_observations_list))
  } else {
    return(data.frame()) # Restituisce un dataframe vuoto se nessuna coppia è stata generata
  }
}


#' Stima i parametri usando il metodo Pairwise Likelihood
estimate_pairwise_likelihood <- function(paired_df, initial_psi = 0.05, max_iter = 150, tolerance = 1e-7) {
  # Restituisci tipi specifici per coerenza nei risultati
  default_return <- list(beta_hat = NA_real_, sigma_sq_hat = NA_real_, psi_hat = NA_real_, 
                         converged = FALSE, iterations = 0L, q_pairs = 0L)
  
  if (is.null(paired_df) || nrow(paired_df) == 0) {
    return(default_return)
  }
  
  q_pairs <- nrow(paired_df)
  if (q_pairs == 0) {
    default_return$q_pairs <- 0L
    return(default_return)
  }
  
  all_x_observed <- c(paired_df$x_i, paired_df$x_l)
  all_y_observed <- c(paired_df$y_i, paired_df$y_l)
  
  if(length(all_x_observed[!is.na(all_x_observed)]) == 0) mean_x <- 0 else mean_x <- mean(all_x_observed, na.rm = TRUE)
  if(length(all_y_observed[!is.na(all_y_observed)]) == 0) mean_y <- 0 else mean_y <- mean(all_y_observed, na.rm = TRUE)
  
  xi_c <- paired_df$x_i - mean_x
  xl_c <- paired_df$x_l - mean_x
  yi_c <- paired_df$y_i - mean_y
  yl_c <- paired_df$y_l - mean_y
  
  alpha_1 <- sum(xi_c^2 + xl_c^2, na.rm = TRUE)
  alpha_2 <- sum(yi_c^2 + yl_c^2, na.rm = TRUE)
  alpha_3 <- sum(xi_c * yi_c + xl_c * yl_c, na.rm = TRUE)
  alpha_4 <- sum(xi_c * yl_c + xl_c * yi_c, na.rm = TRUE)
  alpha_5 <- sum(xi_c * xl_c, na.rm = TRUE)
  alpha_6 <- sum(yi_c * yl_c, na.rm = TRUE)
  
  if(any(is.na(c(alpha_1, alpha_2, alpha_3, alpha_4, alpha_5, alpha_6)))){
    default_return$q_pairs <- q_pairs
    return(default_return)
  }
  
  psi_hat <- initial_psi
  beta_hat <- NA_real_
  sigma_sq_hat <- NA_real_
  
  converged <- FALSE
  iterations_done <- 0L # Tipo Integer
  param_hist <- matrix(NA_real_, nrow = max_iter + 1, ncol = 3) 
  colnames(param_hist) <- c("beta", "sigma_sq", "psi")
  param_hist[1,] <- c(beta_hat, sigma_sq_hat, psi_hat)
  
  for (iter in 1:max_iter) {
    iterations_done <- iter
    beta_hat_old_iter <- param_hist[iter, "beta"]
    sigma_sq_hat_old_iter <- param_hist[iter, "sigma_sq"]
    psi_hat_old_iter <- param_hist[iter, "psi"]
    
    current_psi_for_beta_calc <- ifelse(is.na(psi_hat_old_iter), initial_psi, psi_hat_old_iter)
    denominator_beta <- alpha_1 - 2 * current_psi_for_beta_calc * alpha_5
    if (is.na(denominator_beta) || abs(denominator_beta) < 1e-12) { 
      beta_hat <- ifelse(is.na(beta_hat_old_iter), 0, beta_hat_old_iter) 
    } else {
      beta_hat <- (alpha_3 - current_psi_for_beta_calc * alpha_4) / denominator_beta
    }
    if(is.na(beta_hat)) beta_hat <- ifelse(is.na(beta_hat_old_iter), 0, beta_hat_old_iter) 
    
    current_psi_for_sigma_calc <- ifelse(is.na(psi_hat_old_iter), initial_psi, psi_hat_old_iter)
    numerator_sigma_sq <- alpha_2 + beta_hat^2 * alpha_1 - 2 * beta_hat * alpha_3 -
      2 * current_psi_for_sigma_calc * alpha_6 - 2 * current_psi_for_sigma_calc * beta_hat^2 * alpha_5 +
      2 * current_psi_for_sigma_calc * beta_hat * alpha_4
    denominator_sigma_sq <- 2 * q_pairs * (1 - current_psi_for_sigma_calc^2)
    
    if (is.na(numerator_sigma_sq) || is.na(denominator_sigma_sq) || abs(denominator_sigma_sq) < 1e-12 || numerator_sigma_sq <= 1e-9) { 
      sigma_sq_hat <- ifelse(is.na(sigma_sq_hat_old_iter) || (!is.na(sigma_sq_hat_old_iter) && sigma_sq_hat_old_iter <= 1e-9), 1e-6, sigma_sq_hat_old_iter)
    } else {
      sigma_sq_hat <- numerator_sigma_sq / denominator_sigma_sq
      if (is.na(sigma_sq_hat) || sigma_sq_hat <= 1e-9) { sigma_sq_hat <- 1e-6 } 
    }
    if(is.na(sigma_sq_hat)) sigma_sq_hat <- ifelse(is.na(sigma_sq_hat_old_iter) || sigma_sq_hat_old_iter <= 1e-9, 1e-6, sigma_sq_hat_old_iter)
    
    numerator_psi <- alpha_6 - beta_hat * alpha_4 + beta_hat^2 * alpha_5
    denominator_psi <- q_pairs * sigma_sq_hat 
    
    psi_hat_calc_new <- psi_hat_old_iter 
    if (!is.na(numerator_psi) && !is.na(denominator_psi) && abs(denominator_psi) > 1e-12) { 
      psi_hat_calc_new <- numerator_psi / denominator_psi
    }
    if(is.na(psi_hat_calc_new)) psi_hat_calc_new <- ifelse(is.na(psi_hat_old_iter), initial_psi, psi_hat_old_iter)
    
    psi_hat <- max(-0.999, min(0.999, psi_hat_calc_new)) 
    
    param_hist[iter+1,] <- c(beta_hat, sigma_sq_hat, psi_hat)
    
    if (iter > 1 &&
        !is.na(psi_hat_old_iter) && !is.na(beta_hat_old_iter) && !is.na(sigma_sq_hat_old_iter) &&
        !is.na(psi_hat) && !is.na(beta_hat) && !is.na(sigma_sq_hat)) {
      diff_beta <- abs(beta_hat - beta_hat_old_iter)
      diff_sigma <- abs(sigma_sq_hat - sigma_sq_hat_old_iter)
      diff_psi <- abs(psi_hat - psi_hat_old_iter)
      if (max(c(diff_beta, diff_sigma, diff_psi), na.rm = TRUE) < tolerance) {
        converged <- TRUE
        break
      }
    }
  }
  final_params_idx <- iterations_done + 1 
  return(list(beta_hat = param_hist[final_params_idx, "beta"],
              sigma_sq_hat = param_hist[final_params_idx, "sigma_sq"],
              psi_hat = param_hist[final_params_idx, "psi"],
              converged = converged, iterations = iterations_done, q_pairs = q_pairs))
}


#' Funzione per visualizzare e salvare uno snapshot della simulazione
plot_simulation_snapshot <- function(base_data_sf,
                                     all_candidate_h3_polys_sf,
                                     selected_h3_polys_sf,
                                     paired_observations_df,
                                     title_suffix = "",
                                     save_to_disk = FALSE,
                                     output_dir = "plots",
                                     filename_base = "snapshot") {
  
  cat(paste0("DEBUG plot_simulation_snapshot: Entrato. save_to_disk = ", save_to_disk, 
             ", output_dir = ", output_dir, ", filename_base = ", filename_base, "\n"))
  cat("DEBUG plot_simulation_snapshot: class(base_data_sf): ", class(base_data_sf)[1], ", nrow: ", if(is.null(base_data_sf)) 0 else nrow(base_data_sf), "\n")
  if(!is.null(all_candidate_h3_polys_sf)) cat("DEBUG plot_simulation_snapshot: class(all_candidate_h3_polys_sf): ", class(all_candidate_h3_polys_sf)[1], ", nrow: ", nrow(all_candidate_h3_polys_sf), "\n") else cat("DEBUG plot_simulation_snapshot: all_candidate_h3_polys_sf is NULL\n")
  if(!is.null(selected_h3_polys_sf)) cat("DEBUG plot_simulation_snapshot: class(selected_h3_polys_sf): ", class(selected_h3_polys_sf)[1], ", nrow: ", nrow(selected_h3_polys_sf), "\n") else cat("DEBUG plot_simulation_snapshot: selected_h3_polys_sf is NULL\n")
  
  plot_title <- paste("Snapshot Simulazione:", title_suffix)
  paired_points_sf_from_df <- st_sf(geometry = st_sfc(crs = st_crs(base_data_sf))) # Inizializza vuoto con CRS
  
  if (!is.null(paired_observations_df) && nrow(paired_observations_df) > 0 &&
      "unique_id_i" %in% names(paired_observations_df) &&
      "unique_id_l" %in% names(paired_observations_df)) {
    paired_points_ids <- unique(c(paired_observations_df$unique_id_i, paired_observations_df$unique_id_l))
    if (length(paired_points_ids) > 0) {
      if("unique_id" %in% names(base_data_sf) && nrow(base_data_sf) > 0){ # Aggiunto controllo nrow(base_data_sf)
        paired_points_sf_from_df <- base_data_sf %>% filter(unique_id %in% paired_points_ids)
      } else {
        warning("Colonna 'unique_id' non trovata in base_data_sf o base_data_sf è vuoto per plot_simulation_snapshot.")
      }
    }
  }
  if(!is.null(paired_points_sf_from_df)) cat("DEBUG plot_simulation_snapshot: class(paired_points_sf_from_df): ", class(paired_points_sf_from_df)[1], ", nrow: ", nrow(paired_points_sf_from_df), "\n") else cat("DEBUG plot_simulation_snapshot: paired_points_sf_from_df is NULL (o non creato)\n")
  
  p <- ggplot() 
  
  # Aggiungi layer solo se i dati esistono e sono validi
  if (!is.null(base_data_sf) && nrow(base_data_sf) > 0 && inherits(st_geometry(base_data_sf), "sfc_POINT")) {
    p <- p + geom_sf(data = base_data_sf, aes(color = "Tutti i Punti"), size = 0.3, alpha = 0.2, show.legend = "point")
  } else {
    cat("DEBUG plot_simulation_snapshot: Skipping base_data_sf layer (vuoto o non point).\n")
  }
  
  if (!is.null(all_candidate_h3_polys_sf) && nrow(all_candidate_h3_polys_sf) > 0 && inherits(st_geometry(all_candidate_h3_polys_sf), "sfc_POLYGON")) {
    p <- p + geom_sf(data = all_candidate_h3_polys_sf, aes(fill = "Celle H3 Candidate Iniziali"), alpha = 0.05, color = "grey80", show.legend = "fill")
  } else {
    cat("DEBUG plot_simulation_snapshot: Skipping all_candidate_h3_polys_sf layer.\n")
  }
  
  if (!is.null(selected_h3_polys_sf) && nrow(selected_h3_polys_sf) > 0 && inherits(st_geometry(selected_h3_polys_sf), "sfc_POLYGON")) {
    p <- p + geom_sf(data = selected_h3_polys_sf, aes(fill = "Celle H3 Selezionate Finali"), alpha = 0.2, color = "blue", show.legend = "fill")
  } else {
    cat("DEBUG plot_simulation_snapshot: Skipping selected_h3_polys_sf layer.\n")
  }
  
  if (!is.null(paired_points_sf_from_df) && nrow(paired_points_sf_from_df) > 0 && inherits(st_geometry(paired_points_sf_from_df), "sfc_POINT")) {
    p <- p + geom_sf(data = paired_points_sf_from_df, aes(color = "Punti nelle Coppie"), size = 1.0, shape=21, fill="red", alpha=0.7, show.legend = "point")
  } else {
    cat("DEBUG plot_simulation_snapshot: Skipping paired_points_sf_from_df layer.\n")
  }
  
  p <- p +
    scale_color_manual(name = "Legenda Punti",
                       values = c("Tutti i Punti" = "grey60", "Punti nelle Coppie" = "red"),
                       breaks = c("Tutti i Punti", "Punti nelle Coppie"),
                       drop = FALSE) + # drop = FALSE per mantenere tutte le legende
    scale_fill_manual(name = "Legenda Celle H3",
                      values = c("Celle H3 Candidate Iniziali" = "lightyellow", "Celle H3 Selezionate Finali" = "skyblue"),
                      breaks = c("Celle H3 Candidate Iniziali", "Celle H3 Selezionate Finali"),
                      guide = guide_legend(override.aes = list(alpha = 0.3)),
                      drop = FALSE) + # drop = FALSE
    labs(title = plot_title, x = "Longitudine", y = "Latitudine") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size=9, hjust = 0.5),
          legend.box = "vertical",
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
          legend.spacing.y = unit(0.1, 'cm'))
  
  plot_printed_successfully <- FALSE
  tryCatch({
    # Rimuovi suppressWarnings temporaneamente per vedere tutti i messaggi da ggplot
    print(p) 
    plot_printed_successfully <- TRUE
    cat("DEBUG plot_simulation_snapshot: print(p) eseguito con successo.\n")
  }, error = function(e) {
    cat("Errore durante la generazione del plot snapshot (stampa a schermo):", conditionMessage(e), "\n")
    # print(e) # Per un traceback completo dell'errore
  })
  
  if (save_to_disk && plot_printed_successfully) {
    cat(paste0("DEBUG plot_simulation_snapshot: Condizioni per salvare soddisfatte. Tentativo di salvataggio.\n"))
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      cat(paste0("DEBUG plot_simulation_snapshot: Creata directory: ", output_dir, "\n"))
    }
    plot_filepath <- file.path(output_dir, paste0(filename_base, ".pdf")) 
    cat(paste0("    DEBUG: Tentativo di salvataggio snapshot in: ", plot_filepath, "\n"))
    tryCatch({
      # Usare "pdf" come device standard. "cairo_pdf" è un'opzione se Cairo è configurato.
      ggsave(filename = plot_filepath, plot = p, width = 10, height = 8, device = "pdf", bg = "white") 
      cat(paste0("    Plot snapshot salvato in: ", plot_filepath, "\n"))
    }, error = function(e_save) {
      cat(paste0("    Errore durante il salvataggio del plot snapshot su disco: ", conditionMessage(e_save), "\n"))
      # print(e_save) # Per un traceback completo dell'errore di salvataggio
    })
  } else {
    cat(paste0("DEBUG plot_simulation_snapshot: Condizioni per salvare NON soddisfatte o plot non stampato. save_to_disk=", save_to_disk, 
               ", plot_printed_successfully=", plot_printed_successfully, "\n"))
    if(save_to_disk && !plot_printed_successfully){
      cat(paste0("    ATTENZIONE: Plot non stampato a schermo, quindi non salvato su disco: ", filename_base, ".pdf\n"))
    }
  }
}

#-----------------------------------------------------------------------------#
# 3. DEFINIZIONE DEGLI SCENARI DI SIMULAZIONE ----
#-----------------------------------------------------------------------------#
scenario_parameters <- expand.grid(
  SCENARIO_ID = NA_integer_, 
  TOTAL_POINTS_TO_GENERATE_scenario = c(10000, 25000), 
  H3_RESOLUTION_scenario = c(7),  
  K_RING_SEPARATION_scenario = c(1), 
  PSI_TRUE_scenario = c(0.3), 
  TARGET_Q_PAIRS_scenario = c(50, 100, 250), 
  MIN_OBS_PER_H3_scenario = c(2), 
  DATA_TYPE_scenario = c("uniform", "clustered"),
  BETA_TRUE_scenario = c(1.5),
  SIGMA_SQ_TRUE_scenario = c(1.0), 
  MAX_ITERATIONS_H3_SELECT_scenario = c(200),
  N_SIMULATIONS_per_scenario = 10, 
  K_FOR_GM_WEIGHTS_scenario = c(4, 6), 
  RUN_GM_MODEL_scenario = c(TRUE), 
  stringsAsFactors = FALSE
)
scenario_parameters$SCENARIO_ID <- 1:nrow(scenario_parameters)

X_MEAN_global <- 0
X_SD_global <- 1
IF_CLUSTERED_LAMBDA_CENTROIDS_global <- 5
IF_CLUSTERED_LAMBDA_POINTS_PER_CLUSTER_global <- 100
IF_CLUSTERED_SIGMA_CLUSTER_global <- 0.05

INITIAL_PSI_GUESS_global <- 0.05
NUM_ESTIMATION_ITER_global <- 1000 
CONVERGENCE_TOLERANCE_global <- 1e-7

SET_SEED_global <- 12345
PLOT_FIRST_RUN_SNAPSHOT_global <- TRUE
DEBUG_H3_SELECTION_global <- TRUE 
DEBUG_SCENARIO_ID_H3_SELECTION <- 1:nrow(scenario_parameters) 


#-----------------------------------------------------------------------------#
# 3.5 CONFIG_TAG AUTOMATICO E IMPOSTAZIONI DI SALVATAGGIO ----
#-----------------------------------------------------------------------------#
if (exists("scenario_parameters") && is.data.frame(scenario_parameters)) {
  params_for_hash <- scenario_parameters[, !names(scenario_parameters) %in% c("SCENARIO_ID", "N_SIMULATIONS_per_scenario")]
  CONFIG_TAG <- paste0("ConfigDigest-", digest(params_for_hash, algo = "crc32"))
} else {
  CONFIG_TAG <- "UndefinedConfig_ScenarioParamsNonTrovato"
}

OUTPUT_TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
MAIN_OUTPUT_DIR <- file.path("output", paste0(CONFIG_TAG, "_run_", OUTPUT_TIMESTAMP))
DATA_SAVE_DIR <- file.path(MAIN_OUTPUT_DIR, "data")
PLOT_SNAPSHOT_SAVE_DIR <- file.path(MAIN_OUTPUT_DIR, "plots_snapshot_per_scenario")
PLOT_BOXPLOT_SAVE_DIR <- file.path(MAIN_OUTPUT_DIR, "plots_summary_boxplots")

if(!dir.exists(DATA_SAVE_DIR)) dir.create(DATA_SAVE_DIR, showWarnings = FALSE, recursive = TRUE)
if(!dir.exists(PLOT_SNAPSHOT_SAVE_DIR)) dir.create(PLOT_SNAPSHOT_SAVE_DIR, showWarnings = FALSE, recursive = TRUE)
if(!dir.exists(PLOT_BOXPLOT_SAVE_DIR)) dir.create(PLOT_BOXPLOT_SAVE_DIR, showWarnings = FALSE, recursive = TRUE)

cat(paste0("!!! Output di questa esecuzione (CONFIG_TAG: ", CONFIG_TAG, ") verranno salvati in: ", MAIN_OUTPUT_DIR, " !!!\n\n"))

#-----------------------------------------------------------------------------#
# 4. FUNZIONE PER ESEGUIRE UNO SCENARIO DI SIMULAZIONE ----
#-----------------------------------------------------------------------------#
run_simulation_scenario <- function(params_scenario, sim_extent_polygon) {
  results_list_scenario <- list()
  cat(paste0("    Inizio Scenario ID: ", params_scenario$SCENARIO_ID,
             ", N_Sims: ", params_scenario$N_SIMULATIONS_per_scenario, "\n"))
  
  pb_sims <- progress::progress_bar$new(
    format = paste0("        Simulazione Scenario ", params_scenario$SCENARIO_ID, " [:bar] :percent Eseguite: :current/:total, ETA: :eta"),
    total = params_scenario$N_SIMULATIONS_per_scenario,
    width = 80,
    clear = FALSE
  )
  
  for (sim_run in 1:params_scenario$N_SIMULATIONS_per_scenario) {
    set.seed(SET_SEED_global + params_scenario$SCENARIO_ID * 1000 + sim_run)
    pb_sims$tick()
    
    current_run_results <- data.frame(
      sim_id = sim_run, scenario_id = as.integer(params_scenario$SCENARIO_ID), 
      beta_true = params_scenario$BETA_TRUE_scenario,
      sigma_sq_true = params_scenario$SIGMA_SQ_TRUE_scenario,
      psi_true = params_scenario$PSI_TRUE_scenario,
      beta_hat_pw = NA_real_, sigma_sq_hat_pw = NA_real_, psi_hat_pw = NA_real_,
      converged_pw = FALSE, iterations_pw = NA_integer_, q_pairs_used_pw = 0L,
      time_pw = NA_real_,
      h3_selected_count = 0L,
      beta_hat_gm = NA_real_, lambda_hat_gm = NA_real_, sigma_sq_hat_gm = NA_real_,
      converged_gm = FALSE, time_gm = NA_real_, n_obs_gm = NA_integer_
    )
    
    current_base_data <- generate_base_spatial_data(
      n_points = params_scenario$TOTAL_POINTS_TO_GENERATE_scenario,
      data_type = params_scenario$DATA_TYPE_scenario,
      extent_polygon = sim_extent_polygon,
      x_mean = X_MEAN_global, x_sd = X_SD_global,
      cluster_lambda_centroids = IF_CLUSTERED_LAMBDA_CENTROIDS_global,
      cluster_lambda_points_per_cluster = IF_CLUSTERED_LAMBDA_POINTS_PER_CLUSTER_global,
      cluster_sigma = IF_CLUSTERED_SIGMA_CLUSTER_global
    )
    
    if(nrow(current_base_data) < max(2, params_scenario$MIN_OBS_PER_H3_scenario * 2)) { 
      results_list_scenario[[sim_run]] <- current_run_results
      next
    }
    
    activate_debug_h3 <- DEBUG_H3_SELECTION_global && sim_run == 1 && 
      (length(DEBUG_SCENARIO_ID_H3_SELECTION) == 0 || params_scenario$SCENARIO_ID %in% DEBUG_SCENARIO_ID_H3_SELECTION)
    
    
    all_candidate_polys_sf_for_plot <- NULL
    selected_polys_sf_for_plot <- NULL
    if(PLOT_FIRST_RUN_SNAPSHOT_global && sim_run == 1) {
      temp_base_data_with_h3_for_plot <- current_base_data %>%
        mutate(h3_index_plot = h3jsr::point_to_cell(geometry, res = params_scenario$H3_RESOLUTION_scenario))
      h3_counts_for_plot <- temp_base_data_with_h3_for_plot %>%
        st_drop_geometry() %>% group_by(h3_index_plot) %>%
        summarise(n_obs = n(), .groups = 'drop') %>%
        filter(n_obs >= params_scenario$MIN_OBS_PER_H3_scenario)
      if(nrow(h3_counts_for_plot) > 0 && length(h3_counts_for_plot$h3_index_plot) > 0) {
        all_candidate_polys_sf_for_plot <- cell_to_polygon(h3_counts_for_plot$h3_index_plot, simple = FALSE)
      } else {
        all_candidate_polys_sf_for_plot <- st_sf(geometry = st_sfc(crs = st_crs(current_base_data))) 
      }
      selected_polys_sf_for_plot <- st_sf(geometry = st_sfc(crs = st_crs(current_base_data))) 
    }
    
    current_selected_h3s <- select_h3_cells(
      sf_data = current_base_data,
      h3_res = params_scenario$H3_RESOLUTION_scenario,
      min_obs_in_h3 = params_scenario$MIN_OBS_PER_H3_scenario,
      k_ring_separation = params_scenario$K_RING_SEPARATION_scenario,
      target_selected_h3_count = params_scenario$TARGET_Q_PAIRS_scenario,
      max_iterations_limit = params_scenario$MAX_ITERATIONS_H3_SELECT_scenario,
      scenario_id_for_debug = paste0(params_scenario$SCENARIO_ID, "_Run", sim_run),
      debug_h3_selection = activate_debug_h3
    )
    current_run_results$h3_selected_count <- length(current_selected_h3s)
    
    if (PLOT_FIRST_RUN_SNAPSHOT_global && sim_run == 1 && current_run_results$h3_selected_count > 0) {
      selected_polys_sf_for_plot <- cell_to_polygon(current_selected_h3s, simple = FALSE)
    }
    
    if (current_run_results$h3_selected_count == 0) {
      results_list_scenario[[sim_run]] <- current_run_results
      if (PLOT_FIRST_RUN_SNAPSHOT_global && sim_run == 1) { 
        plot_simulation_snapshot( 
          base_data_sf = current_base_data,
          all_candidate_h3_polys_sf = all_candidate_polys_sf_for_plot,
          selected_h3_polys_sf = selected_polys_sf_for_plot, 
          paired_observations_df = data.frame(), 
          title_suffix = paste0("Scen ", params_scenario$SCENARIO_ID, ", Run ", sim_run, " - NO H3 SELECTED"),
          save_to_disk = TRUE, output_dir = PLOT_SNAPSHOT_SAVE_DIR,
          filename_base = paste0("snapshot_ScenID", params_scenario$SCENARIO_ID,"_Run", sim_run, "_NO_H3_SELECTED")
        )
      }
      next
    }
    
    current_paired_observations_df <- generate_paired_observations(
      all_sf_data = current_base_data,
      selected_h3_ids = current_selected_h3s,
      h3_res = params_scenario$H3_RESOLUTION_scenario,
      beta_true = params_scenario$BETA_TRUE_scenario,
      sigma_sq_true = params_scenario$SIGMA_SQ_TRUE_scenario,
      psi_true = params_scenario$PSI_TRUE_scenario
    )
    
    if (is.null(current_paired_observations_df) || nrow(current_paired_observations_df) == 0) {
      results_list_scenario[[sim_run]] <- current_run_results
      if (PLOT_FIRST_RUN_SNAPSHOT_global && sim_run == 1) { 
        plot_simulation_snapshot( 
          base_data_sf = current_base_data,
          all_candidate_h3_polys_sf = all_candidate_polys_sf_for_plot,
          selected_h3_polys_sf = selected_polys_sf_for_plot, 
          paired_observations_df = data.frame(), 
          title_suffix = paste0("Scen ", params_scenario$SCENARIO_ID, ", Run ", sim_run, " - NO PAIRS"),
          save_to_disk = TRUE, output_dir = PLOT_SNAPSHOT_SAVE_DIR,
          filename_base = paste0("snapshot_ScenID", params_scenario$SCENARIO_ID,"_Run", sim_run, "_NO_PAIRS")
        )
      }
      next
    }
    current_run_results$q_pairs_used_pw <- nrow(current_paired_observations_df)
    
    time_start_pw <- Sys.time()
    estimation_output_pw <- estimate_pairwise_likelihood(
      paired_df = current_paired_observations_df,
      initial_psi = INITIAL_PSI_GUESS_global,
      max_iter = NUM_ESTIMATION_ITER_global,
      tolerance = CONVERGENCE_TOLERANCE_global
    )
    time_end_pw <- Sys.time()
    current_run_results$time_pw <- as.numeric(difftime(time_end_pw, time_start_pw, units = "secs"))
    current_run_results$beta_hat_pw <- estimation_output_pw$beta_hat
    current_run_results$sigma_sq_hat_pw <- estimation_output_pw$sigma_sq_hat
    current_run_results$psi_hat_pw <- estimation_output_pw$psi_hat
    current_run_results$converged_pw <- estimation_output_pw$converged
    current_run_results$iterations_pw <- estimation_output_pw$iterations
    
    if (params_scenario$RUN_GM_MODEL_scenario) {
      ids_i <- current_paired_observations_df$unique_id_i
      ids_l <- current_paired_observations_df$unique_id_l
      all_ids_in_pairs <- unique(c(ids_i, ids_l))
      
      points_for_gm_sf_geom <- current_base_data %>% 
        filter(unique_id %in% all_ids_in_pairs) %>%
        select(unique_id, geometry) 
      
      df_for_gm_data <- data.frame(unique_id = all_ids_in_pairs) 
      
      map_i <- current_paired_observations_df %>% select(unique_id = unique_id_i, y = y_i, x_value = x_i)
      map_l <- current_paired_observations_df %>% select(unique_id = unique_id_l, y = y_l, x_value = x_l)
      map_all_yx <- bind_rows(map_i, map_l) %>% distinct(unique_id, .keep_all = TRUE)
      
      df_for_gm_data <- df_for_gm_data %>%
        left_join(map_all_yx, by = "unique_id")
      
      df_for_gm_sf <- points_for_gm_sf_geom %>%
        left_join(df_for_gm_data, by = "unique_id") %>%
        filter(!is.na(y) & !is.na(x_value) & sf::st_is_valid(geometry)) 
      
      current_run_results$n_obs_gm <- nrow(df_for_gm_sf)
      
      actual_k_for_gm <- min(params_scenario$K_FOR_GM_WEIGHTS_scenario, nrow(df_for_gm_sf) - 1)
      
      if (nrow(df_for_gm_sf) >= actual_k_for_gm + 1 && actual_k_for_gm > 0 && nrow(df_for_gm_sf) > 2 && ncol(st_coordinates(df_for_gm_sf)) == 2) {
        time_start_gm <- Sys.time()
        gm_model <- NULL
        tryCatch({
          coords_gm <- st_coordinates(df_for_gm_sf)
          knn <- knearneigh(coords_gm, k = actual_k_for_gm, longlat = FALSE) 
          nb_gm <- knn2nb(knn, sym = TRUE) # Aggiunto sym=TRUE per sicurezza, anche se k-nn non è simmetrico di per sè
          
          has_neighbours <- card(nb_gm) > 0
          
          if(sum(has_neighbours) <= actual_k_for_gm +1){ 
            current_run_results$converged_gm <- FALSE
          } else {
            # Subset per nodi con vicini per nb2listw e per il modello
            df_for_gm_sf_filtered <- df_for_gm_sf[has_neighbours, ] 
            nb_gm_filtered <- subset(nb_gm, has_neighbours)
            listw_gm <- nb2listw(nb_gm_filtered, style="W", zero.policy = TRUE) 
            
            if(nrow(df_for_gm_sf_filtered) > actual_k_for_gm + 1) { # Ricontrolla dopo il filtro
              gm_model <- GMerrorsar(y ~ x_value, data = df_for_gm_sf_filtered, listw = listw_gm, zero.policy = TRUE, method="nlminb")
              
              coefs_gm <- coef(gm_model)
              current_run_results$beta_hat_gm <- if("x_value" %in% names(coefs_gm)) coefs_gm["x_value"] else NA_real_
              current_run_results$lambda_hat_gm <- gm_model$lambda
              current_run_results$sigma_sq_hat_gm <- gm_model$s2
              current_run_results$converged_gm <- TRUE # GMerrorsar non ha un flag di convergenza esplicito, assumiamo se non dà errore
            } else {
              current_run_results$converged_gm <- FALSE
            }
          }
        }, error = function(e) {
          current_run_results$converged_gm <- FALSE
        })
        time_end_gm <- Sys.time()
        current_run_results$time_gm <- if (current_run_results$converged_gm && !is.null(gm_model)) as.numeric(difftime(time_end_gm, time_start_gm, units = "secs")) else NA_real_
      } else {
        current_run_results$converged_gm <- FALSE
      }
    }
    
    results_list_scenario[[sim_run]] <- current_run_results
    
    if (PLOT_FIRST_RUN_SNAPSHOT_global && sim_run == 1) {
      # ... (codice snapshot come prima)
      plot_title_sfx <- paste0("Scen ", params_scenario$SCENARIO_ID,
                               ", Run ", sim_run,
                               "\nH3Res:", params_scenario$H3_RESOLUTION_scenario,
                               ", KRing:", params_scenario$K_RING_SEPARATION_scenario,
                               ", TargetH3Cells:", params_scenario$TARGET_Q_PAIRS_scenario,
                               "\nActualH3Selected: ", current_run_results$h3_selected_count,
                               ", ActualQpairsPW: ", current_run_results$q_pairs_used_pw,
                               "\nPsiTrue: ", params_scenario$PSI_TRUE_scenario,
                               ", PsiHatPW: ", round(current_run_results$psi_hat_pw,3))
      
      snapshot_filename_base <- paste0(
        "snapshot_ScenID", params_scenario$SCENARIO_ID,
        "_Run", sim_run,
        "_H3Res", params_scenario$H3_RESOLUTION_scenario,
        "_KRing", params_scenario$K_RING_SEPARATION_scenario,
        "_TargetH3Cells", params_scenario$TARGET_Q_PAIRS_scenario,
        "_Npts", params_scenario$TOTAL_POINTS_TO_GENERATE_scenario
      )
      plot_simulation_snapshot(
        base_data_sf = current_base_data,
        all_candidate_h3_polys_sf = all_candidate_polys_sf_for_plot,
        selected_h3_polys_sf = selected_polys_sf_for_plot,
        paired_observations_df = current_paired_observations_df,
        title_suffix = plot_title_sfx,
        save_to_disk = TRUE,
        output_dir = PLOT_SNAPSHOT_SAVE_DIR,
        filename_base = snapshot_filename_base
      )
    }
  } 
  
  cat(paste0("    Scenario ID: ", params_scenario$SCENARIO_ID, " completato.\n"))
  if (length(results_list_scenario) > 0) {
    return(do.call(rbind, results_list_scenario))
  } else {
    return(data.frame(sim_id=integer(), scenario_id=integer(), beta_true=numeric(), sigma_sq_true=numeric(), psi_true=numeric(),
                      beta_hat_pw=NA_real_, sigma_sq_hat_pw=NA_real_, psi_hat_pw=NA_real_,
                      converged_pw=FALSE, iterations_pw=NA_integer_, q_pairs_used_pw=0L, time_pw=NA_real_,
                      h3_selected_count=0L,
                      beta_hat_gm=NA_real_, lambda_hat_gm=NA_real_, sigma_sq_hat_gm=NA_real_,
                      converged_gm=FALSE, time_gm=NA_real_, n_obs_gm=NA_integer_
    ))
  }
}


#-----------------------------------------------------------------------------#
# 5. CICLO PRINCIPALE SU TUTTI GLI SCENARI ----
#-----------------------------------------------------------------------------#
all_scenario_results_list <- list()
coords_extent_main <- matrix(c(0,0, 1,0, 1,1, 0,1, 0,0), ncol=2, byrow=TRUE)
sim_extent_polygon_main <- st_polygon(list(coords_extent_main)) %>% st_sfc(crs = 4326) 

cat(paste0("\nInizio esecuzione di ", nrow(scenario_parameters), " scenari...\n"))

for (i in 1:nrow(scenario_parameters)) {
  current_scenario_params_row <- scenario_parameters[i, ]
  scenario_results_for_this_config_df <- run_simulation_scenario(current_scenario_params_row, sim_extent_polygon_main)
  all_scenario_results_list[[i]] <- scenario_results_for_this_config_df
}
cat("\nEsecuzione di tutti gli scenari completata.\n") # CORRETTO L'ERRORE DI PARENTESI


if (length(all_scenario_results_list) > 0) {
  final_results_df <- tryCatch({
    bind_rows(all_scenario_results_list) 
  }, error = function(e) {
    cat("Errore in bind_rows(all_scenario_results_list): ", e$message, "\n")
    cat("Ispezionare all_scenario_results_list per dataframe con strutture diverse.\n")
    return(data.frame())
  })
} else {
  final_results_df <- data.frame()
}

cat("\n--- DIAGNOSI: Nomi colonne in final_results_df ---\n")
if(nrow(final_results_df) > 0) print(names(final_results_df)) else cat("final_results_df è vuoto.\n")
cat("\n--- DIAGNOSI: Struttura di final_results_df (prime righe) ---\n")
if(nrow(final_results_df) > 0) print(head(final_results_df)) else cat("final_results_df è vuoto.\n")


if (nrow(final_results_df) > 0 && !"scenario_id" %in% names(final_results_df)) {
  warning("Colonna 'scenario_id' mancante in final_results_df prima del merge!")
}


if (nrow(final_results_df) > 0) {
  final_results_df$scenario_id <- as.integer(final_results_df$scenario_id)
  scenario_parameters$SCENARIO_ID <- as.integer(scenario_parameters$SCENARIO_ID)
  
  final_results_df_merged <- final_results_df %>%
    left_join(scenario_parameters, by = c("scenario_id" = "SCENARIO_ID"), suffix = c("", ".param_drop")) %>% 
    select(-ends_with(".param_drop")) 
} else {
  final_results_df_merged <- data.frame()
}

cat("\n--- DIAGNOSI: Nomi colonne in final_results_df_merged ---\n")
if(nrow(final_results_df_merged)>0) print(names(final_results_df_merged)) else cat("final_results_df_merged è vuoto.\n")
cat("\n--- DIAGNOSI: Struttura di final_results_df_merged (prime righe) ---\n")
if(nrow(final_results_df_merged)>0) print(head(final_results_df_merged)) else cat("final_results_df_merged è vuoto.\n")

# --- SALVATAGGIO DEI DATI FINALI ---
if (nrow(final_results_df_merged) > 0) {
  results_rds_filename <- file.path(DATA_SAVE_DIR, paste0("final_results_merged_", OUTPUT_TIMESTAMP, ".rds"))
  tryCatch({
    write_rds(final_results_df_merged, results_rds_filename)
    cat(paste0("\nRisultati completi della simulazione salvati in: ", results_rds_filename, "\n"))
  }, error = function(e){cat(paste0("Errore nel salvare i risultati RDS: ", e$message, "\n"))})
} else {
  cat("\nNessun risultato finale da salvare (final_results_df_merged è vuoto).\n")
}

scenarios_rds_filename <- file.path(DATA_SAVE_DIR, paste0("scenario_parameters_definition_", OUTPUT_TIMESTAMP, ".rds"))
tryCatch({
  write_rds(scenario_parameters, scenarios_rds_filename)
  cat(paste0("Definizione dei parametri degli scenari usata salvata in: ", scenarios_rds_filename, "\n"))
}, error = function(e){cat(paste0("Errore nel salvare i parametri scenario RDS: ", e$message, "\n"))})

#-----------------------------------------------------------------------------#
# 6. ANALISI DEI RISULTATI (PER SCENARIO) ----
#-----------------------------------------------------------------------------#
# ... (come prima, con la logica di summarise e plotting aggiornata) ...
if (exists("final_results_df_merged") && !is.null(final_results_df_merged) && nrow(final_results_df_merged) > 0) {
  
  results_df_clean <- final_results_df_merged 
  
  cat("\n--- DIAGNOSI (SEZIONE 6): Nomi colonne in results_df_clean prima di summarise ---\n")
  if(nrow(results_df_clean)>0) print(names(results_df_clean)) else cat("results_df_clean è vuoto.\n")
  cat("\n--- DIAGNOSI (SEZIONE 6): Struttura di results_df_clean (prime righe) prima di summarise ---\n")
  if(nrow(results_df_clean)>0) print(head(results_df_clean)) else cat("results_df_clean è vuoto.\n")
  
  if (nrow(results_df_clean) > 0) {
    param_cols_from_scenario_df <- names(scenario_parameters)[names(scenario_parameters) != "SCENARIO_ID"]
    
    summary_stats_by_scenario <- results_df_clean %>%
      group_by(scenario_id) %>% 
      summarise(
        # Recupera i parametri dello scenario
        across(all_of(param_cols_from_scenario_df), ~first(na.omit(.))),
        
        N_Valid_Runs_PW = sum(!is.na(beta_hat_pw) & !is.na(sigma_sq_hat_pw) & !is.na(psi_hat_pw), na.rm = TRUE),
        N_Converged_PW = sum(converged_pw, na.rm = TRUE),
        Mean_Time_PW = mean(time_pw, na.rm = TRUE),
        
        Mean_Beta_Hat_PW = mean(beta_hat_pw, na.rm = TRUE),
        Median_Beta_Hat_PW = median(beta_hat_pw, na.rm = TRUE),
        SD_Beta_Hat_PW = sd(beta_hat_pw, na.rm = TRUE),
        Bias_Beta_PW = if(any(!is.na(beta_hat_pw))) mean(beta_hat_pw - BETA_TRUE_scenario, na.rm = TRUE) else NA_real_,
        MSE_Beta_PW = if(any(!is.na(beta_hat_pw))) mean((beta_hat_pw - BETA_TRUE_scenario)^2, na.rm = TRUE) else NA_real_,
        
        Mean_SigmaSq_Hat_PW = mean(sigma_sq_hat_pw, na.rm = TRUE),
        Median_SigmaSq_Hat_PW = median(sigma_sq_hat_pw, na.rm = TRUE),
        SD_SigmaSq_Hat_PW = sd(sigma_sq_hat_pw, na.rm = TRUE),
        Bias_SigmaSq_PW = if(any(!is.na(sigma_sq_hat_pw))) mean(sigma_sq_hat_pw - SIGMA_SQ_TRUE_scenario, na.rm = TRUE) else NA_real_,
        MSE_SigmaSq_PW = if(any(!is.na(sigma_sq_hat_pw))) mean((sigma_sq_hat_pw - SIGMA_SQ_TRUE_scenario)^2, na.rm = TRUE) else NA_real_,
        
        Mean_Psi_Hat_PW = mean(psi_hat_pw, na.rm = TRUE),
        Median_Psi_Hat_PW = median(psi_hat_pw, na.rm = TRUE),
        SD_Psi_Hat_PW = sd(psi_hat_pw, na.rm = TRUE),
        Bias_Psi_PW = if(any(!is.na(psi_hat_pw))) mean(psi_hat_pw - PSI_TRUE_scenario, na.rm = TRUE) else NA_real_,
        MSE_Psi_PW = if(any(!is.na(psi_hat_pw))) mean((psi_hat_pw - PSI_TRUE_scenario)^2, na.rm = TRUE) else NA_real_,
        
        Avg_Q_Pairs_Used_PW = mean(q_pairs_used_pw, na.rm = TRUE),
        Avg_H3_Selected_Count = mean(h3_selected_count, na.rm = TRUE),
        ConvergenceRate_PW = if(sum(!is.na(beta_hat_pw) & !is.na(sigma_sq_hat_pw) & !is.na(psi_hat_pw)) > 0) sum(converged_pw, na.rm = TRUE) / sum(!is.na(beta_hat_pw) & !is.na(sigma_sq_hat_pw) & !is.na(psi_hat_pw)) else 0,
        Avg_Iterations_To_Converge_PW = mean(iterations_pw[which(converged_pw)], na.rm = TRUE), 
        
        N_Valid_Runs_GM = if(first(RUN_GM_MODEL_scenario)) sum(!is.na(beta_hat_gm) & !is.na(lambda_hat_gm) & !is.na(sigma_sq_hat_gm), na.rm = TRUE) else 0L,
        N_Converged_GM = if(first(RUN_GM_MODEL_scenario)) sum(converged_gm, na.rm = TRUE) else 0L, 
        Mean_Time_GM = if(first(RUN_GM_MODEL_scenario)) mean(time_gm, na.rm = TRUE) else NA_real_,
        
        Mean_Beta_Hat_GM = if(first(RUN_GM_MODEL_scenario)) mean(beta_hat_gm, na.rm = TRUE) else NA_real_,
        Median_Beta_Hat_GM = if(first(RUN_GM_MODEL_scenario)) median(beta_hat_gm, na.rm = TRUE) else NA_real_,
        SD_Beta_Hat_GM = if(first(RUN_GM_MODEL_scenario)) sd(beta_hat_gm, na.rm = TRUE) else NA_real_,
        Bias_Beta_GM = if(first(RUN_GM_MODEL_scenario) && any(!is.na(beta_hat_gm))) mean(beta_hat_gm - BETA_TRUE_scenario, na.rm = TRUE) else NA_real_,
        MSE_Beta_GM = if(first(RUN_GM_MODEL_scenario) && any(!is.na(beta_hat_gm))) mean((beta_hat_gm - BETA_TRUE_scenario)^2, na.rm = TRUE) else NA_real_,
        
        Mean_Lambda_Hat_GM = if(first(RUN_GM_MODEL_scenario)) mean(lambda_hat_gm, na.rm = TRUE) else NA_real_,
        Median_Lambda_Hat_GM = if(first(RUN_GM_MODEL_scenario)) median(lambda_hat_gm, na.rm = TRUE) else NA_real_,
        SD_Lambda_Hat_GM = if(first(RUN_GM_MODEL_scenario)) sd(lambda_hat_gm, na.rm = TRUE) else NA_real_,
        
        Mean_SigmaSq_Hat_GM = if(first(RUN_GM_MODEL_scenario)) mean(sigma_sq_hat_gm, na.rm = TRUE) else NA_real_,
        Median_SigmaSq_Hat_GM = if(first(RUN_GM_MODEL_scenario)) median(sigma_sq_hat_gm, na.rm = TRUE) else NA_real_,
        SD_SigmaSq_Hat_GM = if(first(RUN_GM_MODEL_scenario)) sd(sigma_sq_hat_gm, na.rm = TRUE) else NA_real_,
        Bias_SigmaSq_GM = if(first(RUN_GM_MODEL_scenario) && any(!is.na(sigma_sq_hat_gm))) mean(sigma_sq_hat_gm - SIGMA_SQ_TRUE_scenario, na.rm = TRUE) else NA_real_, 
        MSE_SigmaSq_GM = if(first(RUN_GM_MODEL_scenario) && any(!is.na(sigma_sq_hat_gm))) mean((sigma_sq_hat_gm - SIGMA_SQ_TRUE_scenario)^2, na.rm = TRUE) else NA_real_,
        
        .groups = 'drop'
      )
    
    cat("\n--- Statistiche Riassuntive per Scenario ---\n")
    print(as.data.frame(summary_stats_by_scenario))
    
    # ... (salvataggio CSV e RDS come prima) ...
    summary_csv_filename <- file.path(DATA_SAVE_DIR, paste0("summary_stats_by_scenario_", OUTPUT_TIMESTAMP, ".csv"))
    tryCatch({
      write_csv(summary_stats_by_scenario, summary_csv_filename)
      cat(paste0("\nStatistiche riassuntive per scenario salvate in formato CSV: ", summary_csv_filename, "\n"))
    }, error = function(e){cat(paste0("Errore nel salvare le statistiche CSV: ", e$message, "\n"))})
    
    summary_rds_filename <- file.path(DATA_SAVE_DIR, paste0("summary_stats_by_scenario_", OUTPUT_TIMESTAMP, ".rds"))
    tryCatch({
      write_rds(summary_stats_by_scenario, summary_rds_filename)
      cat(paste0("Statistiche riassuntive per scenario salvate in formato RDS: ", summary_rds_filename, "\n"))
    }, error = function(e){cat(paste0("Errore nel salvare le statistiche RDS: ", e$message, "\n"))})
    
    
    #-----------------------------------------------------------------------------#
    # 7. VISUALIZZAZIONE (BOXPLOT PER SCENARIO) ----
    #-----------------------------------------------------------------------------#
    scenario_labels_df_for_plot <- scenario_parameters %>%
      mutate(Scenario_Label = paste0("ScenID:", SCENARIO_ID, 
                                     "\nNpts:", TOTAL_POINTS_TO_GENERATE_scenario,
                                     " H3Res:", H3_RESOLUTION_scenario,
                                     " KRing:", K_RING_SEPARATION_scenario,
                                     "\nPsiTrue:", PSI_TRUE_scenario,
                                     " TargetH3:", TARGET_Q_PAIRS_scenario,
                                     " MinObsH3:", MIN_OBS_PER_H3_scenario
      )) %>% select(SCENARIO_ID, Scenario_Label)
    
    
    results_df_for_plot <- results_df_clean %>%
      left_join(scenario_labels_df_for_plot, by = c("scenario_id" = "SCENARIO_ID"))
    
    
    results_long_pw <- results_df_for_plot %>%
      filter(!is.na(beta_hat_pw) & !is.na(sigma_sq_hat_pw) & !is.na(psi_hat_pw)) %>%  # Filtra run validi per PW
      select(scenario_id, Scenario_Label, sim_id, BETA_TRUE_scenario, SIGMA_SQ_TRUE_scenario, PSI_TRUE_scenario, 
             beta_hat_pw, sigma_sq_hat_pw, psi_hat_pw) %>%
      pivot_longer(cols = c(beta_hat_pw, sigma_sq_hat_pw, psi_hat_pw),
                   names_to = "parameter_pw", values_to = "estimate_pw") %>%
      mutate(
        true_value = case_when(
          parameter_pw == "beta_hat_pw" ~ BETA_TRUE_scenario,
          parameter_pw == "sigma_sq_hat_pw" ~ SIGMA_SQ_TRUE_scenario,
          parameter_pw == "psi_hat_pw" ~ PSI_TRUE_scenario 
        ),
        parameter_pw_label = factor(case_when(
          parameter_pw == "beta_hat_pw" ~ "Beta (Pairwise)",
          parameter_pw == "sigma_sq_hat_pw" ~ "Sigma^2 (Pairwise)",
          parameter_pw == "psi_hat_pw" ~ "Psi (Pairwise)"
        ), levels = c("Beta (Pairwise)", "Sigma^2 (Pairwise)", "Psi (Pairwise)"))
      )
    
    if (nrow(results_long_pw) > 0) {
      p_pw_params <- ggplot(results_long_pw, aes(x = Scenario_Label, y = estimate_pw)) +
        geom_boxplot(aes(fill = Scenario_Label), show.legend = FALSE, na.rm = TRUE, outlier.size = 0.5) +
        geom_hline(aes(yintercept = true_value, color = "Valore Vero"), linetype = "dashed") +
        scale_color_manual(name = "", values = c("Valore Vero" = "red")) +
        facet_wrap(~parameter_pw_label, scales = "free_y", ncol = 3) + 
        labs(title = "Distribuzione Stime Parametri Pairwise per Scenario",
             x = "Scenario", y = "Valore Stimato") +
        theme_bw(base_size = 9) +
        theme(axis.text.x = element_text(angle = 65, hjust = 1, size=7),
              legend.position = "top", plot.title = element_text(hjust = 0.5, size=11),
              strip.text = element_text(size=9),
              panel.spacing = unit(0.8, "lines"))
      
      pw_boxplot_filename <- file.path(PLOT_BOXPLOT_SAVE_DIR, paste0("boxplot_pairwise_params_", OUTPUT_TIMESTAMP, ".pdf"))
      tryCatch({
        n_scen_pw <- length(unique(results_long_pw$Scenario_Label))
        plot_width_pw <- max(7, 2.5 * n_scen_pw) 
        plot_height_pw <- max(5, 4 * 1) 
        ggsave(filename = pw_boxplot_filename, plot = p_pw_params,
               width = plot_width_pw, height = plot_height_pw, bg = "white", limitsize = FALSE) # Rimosso cairo_pdf per ora
        cat(paste0("Boxplot parametri Pairwise salvato in: ", pw_boxplot_filename, "\n"))
      }, error = function(e){cat(paste0("Errore nel salvare boxplot Pairwise: ", e$message, "\n"))})
    }
    
    results_long_gm <- results_df_for_plot %>%
      filter(RUN_GM_MODEL_scenario & !is.na(beta_hat_gm) & !is.na(lambda_hat_gm) & !is.na(sigma_sq_hat_gm)) %>% # Filtra run validi per GM
      select(scenario_id, Scenario_Label, sim_id, BETA_TRUE_scenario, SIGMA_SQ_TRUE_scenario, 
             beta_hat_gm, lambda_hat_gm, sigma_sq_hat_gm) %>%
      pivot_longer(cols = c(beta_hat_gm, lambda_hat_gm, sigma_sq_hat_gm),
                   names_to = "parameter_gm", values_to = "estimate_gm") %>%
      mutate(
        true_value_gm = case_when( 
          parameter_gm == "beta_hat_gm" ~ BETA_TRUE_scenario,
          parameter_gm == "sigma_sq_hat_gm" ~ SIGMA_SQ_TRUE_scenario, 
          TRUE ~ NA_real_ 
        ),
        parameter_gm_label = factor(case_when(
          parameter_gm == "beta_hat_gm" ~ "Beta (GM)",
          parameter_gm == "lambda_hat_gm" ~ "Lambda (GM)",
          parameter_gm == "sigma_sq_hat_gm" ~ "Sigma^2 (GM)"
        ), levels = c("Beta (GM)", "Lambda (GM)", "Sigma^2 (GM)"))
      )
    
    if (nrow(results_long_gm) > 0) {
      p_gm_params <- ggplot(results_long_gm, aes(x = Scenario_Label, y = estimate_gm)) +
        geom_boxplot(aes(fill = Scenario_Label), show.legend = FALSE, na.rm = TRUE, outlier.size = 0.5) +
        geom_hline(data = . %>% filter(!is.na(true_value_gm)), aes(yintercept = true_value_gm, color = "Valore Vero Riferimento (PW)"), linetype = "dashed") +
        scale_color_manual(name = "", values = c("Valore Vero Riferimento (PW)" = "red")) +
        facet_wrap(~parameter_gm_label, scales = "free_y", ncol = 3) + 
        labs(title = "Distribuzione Stime Parametri GMerrorsar per Scenario",
             x = "Scenario", y = "Valore Stimato") +
        theme_bw(base_size = 9) +
        theme(axis.text.x = element_text(angle = 65, hjust = 1, size=7),
              legend.position = "top", plot.title = element_text(hjust = 0.5, size=11),
              strip.text = element_text(size=9),
              panel.spacing = unit(0.8, "lines"))
      
      gm_boxplot_filename <- file.path(PLOT_BOXPLOT_SAVE_DIR, paste0("boxplot_gm_params_", OUTPUT_TIMESTAMP, ".pdf"))
      tryCatch({
        n_scen_gm <- length(unique(results_long_gm$Scenario_Label))
        plot_width_gm <- max(7, 2.5 * n_scen_gm) 
        plot_height_gm <- max(5, 4 * 1) 
        ggsave(filename = gm_boxplot_filename, plot = p_gm_params,
               width = plot_width_gm, height = plot_height_gm, bg = "white", limitsize = FALSE) # Rimosso cairo_pdf
        cat(paste0("Boxplot parametri GM salvato in: ", gm_boxplot_filename, "\n"))
      }, error = function(e){cat(paste0("Errore nel salvare boxplot GM: ", e$message, "\n"))})
    }
    
  } else {
    cat("Nessun run di simulazione valido per l'analisi e la visualizzazione.\n")
  }
} else {
  cat("Nessun risultato di simulazione prodotto (final_results_df_merged è vuoto).\n")
}

cat("\n--- Script di simulazione completato ---\n")

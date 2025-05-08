library(sf)
library(dplyr)
library(spatialreg)
library(spdep)
library(h3jsr)


h3_res <- 7

data <- read_sf("kingcounty/kc_house.shp") %>% 
  mutate(h3=point_to_cell(., res=h3_res)) 

h3_to_keep <- data %>%
  group_by(h3) %>% 
  summarise(n=n()) %>%
  filter(n>2) %>% 
  pull(h3)

data_full <- data %>% 
  filter(h3 %in% h3_to_keep)

data_h3 <- cell_to_polygon(h3_to_keep, simple = FALSE) %>%
  st_transform(4326)

ggplot(data=data_h3) +
  geom_sf() +
  geom_sf(data=data)

results <- list()
for(sim in 1:10){
data <- data_full
h3_list <- h3_to_keep
rep <- 0
l <- length(h3_list)
selected <- c()

while(rep < 5){
  s <- sample(h3_list, 1)
  selected <- c(selected, s)
  to_remove <- get_ring(s, ring_size = 1, simple = TRUE)[[1]]
  h3_list <- h3_list[!(h3_list %in% to_remove)]
  
  if(length(h3_list) == l) rep = rep+1 else rep=0
  l = length(h3_list)
  print(rep)
}



data <- data %>% 
  filter(h3 %in% selected) %>% 
  group_by(h3) %>% 
  slice_sample(n=2)


data_h3 %>% 
  filter(h3_address %in% selected) %>% 
  ggplot() +
  geom_sf() +
  geom_sf(data=data)




# soluzione più rapida imponendo block diagonal by 2
grid <- expand.grid(x=1:nrow(data), y=1:nrow(data))
grid <- as.matrix(grid)
match <- apply(grid, 1, function(x) as.numeric(data$h3[x[1]] == data$h3[x[2]]))
m <- matrix(match, nrow=nrow(data))
diag(m) <- 0


X_centered <- data$sqft_liv - mean(data$sqft_liv)
Y_centered <- data$price - mean(data$price)


coppie <- grid[match == 1,]
coppie <- coppie[coppie[,1] != coppie[,2],]
coppie <- coppie[coppie[,1] < coppie[,2],]
coppie
q_pairs <- nrow(coppie)
coppie <- data.frame(i=coppie[,1], l=coppie[,2], 
           xi=X_centered[coppie[,1]], 
           xl=X_centered[coppie[,2]], 
           yi=Y_centered[coppie[,1]], 
           yl=Y_centered[coppie[,2]])

alpha_1 <- sum(coppie$xi^2 + coppie$xl^2)
alpha_2 <- sum(coppie$yi^2 + coppie$yl^2)
alpha_3 <- sum(coppie$xi*coppie$yi + coppie$xl*coppie$yl)
alpha_4 <- sum(coppie$xi*coppie$yl + coppie$xl*coppie$yi)
alpha_5 <- sum(coppie$xi*coppie$xl)
alpha_6 <- sum(coppie$yi*coppie$yl)


verbose = T
psi_hat <- 0
num_estimation_iter <- 100
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

results[[sim]] <- c(beta_hat, sigma_sq_hat, psi_hat)
print(sim)
}

results

mean(sapply(results, function(x) x[1]))
mean(sapply(results, function(x) x[2]))
mean(sapply(results, function(x) x[3]))

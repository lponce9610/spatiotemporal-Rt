# Packages -----------------
library(dplyr)
library(parallel)
library(truncnorm)  # For truncated normal distribution
library(sf)
library(rnaturalearth) # needs rnaturalearthhires for scale = "large"
library(simcausal) # for rbern
library(purrr)
library(stringr)
library(caret) # for sensitivity/specificity
library(pROC)

# Set simulation parameters --------------------
# Functions for temperature-dependent generation time
# Function for parameter in Extrinsic Incubation Period
lambdaEIP<-function(T,v=4.3,beta0=7.9,betat=-0.21,Tbar=0) v/exp(beta0+betat*(T-Tbar))  

# Function for parameter in Intrinsic Incubation Period
lambdaIIP<-function(v=16,beta0=1.78) v/exp(beta0) 

# Computes gamma_k value used in the formula
gammak <- function (a,b,k) {
  b1 <- min(b)
  lga <- log(a)+ k*log(1-b1/b) - log(k)
  return (sum(exp(lga)))
}

# implements the Moschopolous method for finding the probability density
# y: the input value
# a: shape
# b: scale
sum_gamma_dist <- function (y, a, b, K=100){
  b1 <- min(b)
  rho <- sum(a)
  C <- exp(sum(a*log(b1/b)))
  
  delta <- gammav <- c(1,rep(0,K))
  
  gammav <- sapply(1:K, function(x) {gammak(a,b,x)})
  #for (k in 1:K) {
  #  gammav[k] <- gammak(a,b,k)
  #}
  
  delta[1] <- 1
  
  for (k in 0:(K-1)) {
    ss2 <- (1:(k+1))*gammav[1:(k+1)]*delta[1+k+1 - (1:(k+1))]
    delta[(k+1)+1] <- sum(ss2)/(k+1)
  }
  
  xx <- log(delta) + 
    (rho-1 + 1:(K+1))*log(y) -(y/b1) - lgamma(rho+1:(K+1)) - 
    (rho + 1:(K+1))*log(b1)
  gy <- C*sum(exp(xx)) 
  gy
}

# Function to compute the probability distribution, 
# given by \int_0^{y} f(x) dx
# y: the input value
# a: shape
# b: scale
int_sum_gamma <- function(y, a, b, K=200, step=.1, max=50) {
  summ <- 0
  #step <- 0.1
  res <- rep(0, max/step)
  for (i in (1:(max/step))*step) {
    xx <- sum_gamma_dist(i,a,b, K)*step
    if (is.na(xx)) {
      cat("NA i a b K", i,a,b,K)
    }
    if (is.infinite(xx)) {
      cat("infinite i a b K", i,a,b,K)
    }
    if (is.nan(xx)) {
      cat("nan i a b K", i,a,b,K)
    }
    #print(xx)
    summ <- summ + xx
    res[round(i/step)] <- summ
  }
  list(yres=res[round(y/step)], dist=res, i=(1:(max/step))*step)
}

# inverse function given a probability distribution, i.e.,
# given the probability, returns the random variable value
# OBS: it performs a binary search
# p: the probability value
# a: shape
# b: scale
t_sum_gamma_v3 <- function(p, a, b, K=200, step=.1, max=100, res=.1) {  
  xa <- 0
  xb <- max
  while ((xb-xa)>res) {
    y <- (xa+xb)/2.0
    px <- int_sum_gamma(y, a, b, max=max)    
    
    if (px$yres>p) {xb <- y}
    else {xa <- y}
  }
  y
}

# Function to compute the probability distribution, 
# given by \int_0^{y} f(x) dx
# and using time series of temperature values
# y: the input value
# a: shape
# b: scale
# Temp: time series of temperature values
# t: time when Temp[t+1] start
# The idea is to have
# P( X <= t) = P(X > t_i) * P(X > t | X> t_i),
# where t_i is the instant in which temperature has changed previous to 
# time t
# Then, P(X>t | X > t_i) has a time rescaling
# P(X>t | X > t_i) = 1 - \int_{t_{i}}^{t} f(\tau + t_{equiv} - t_{i}, \mathbf{\theta}, \mathbf{\beta(Temp_i)})/P(X > t_i)
# The procedure can continue iteratively
int_sum_gamma_T <- function(y, a, b, Temp, t, K=200, step=.1, max=50, withbreak=TRUE, unitscale=1) {
  summ <- 0
  c <- 1
  b[2] <- 1/lambdaEIP(T=Temp[c])
  tdiff <- 0
  res <- rep(0, max/step)
  pdf <- rep(0, max+1)
  for (i in (1:(max/step))*step) {
    if ((withbreak) && (summ>0.999)) {
      j <- ceiling(i)
      pdf[j+1] <- pdf[j+1] + 0.0
      res[round(i/step)] <- summ        
    }
    else {
      b[2] <- 1/lambdaEIP(T=Temp[c])
      if (c<=length(t)) {
        if (i>t[c]) {
          c <- c+1
          b[2] <- 1/lambdaEIP(T=Temp[c])  
          tsum <- t_sum_gamma_v3(summ, a, b/unitscale, max=max)
          tdiff <- tsum - i
        }
      }    
      xx <- sum_gamma_dist(max(i+tdiff, 0),a,b/unitscale, K)*(step)
      
      summ <- summ + xx
      j <- ceiling(i)
      pdf[j+1] <- pdf[j+1] + xx
      res[round(i/step)] <- summ  
    }
  }
  list(yres=res[y/step], dist=res, i=(1:(max/step))*step, pdf = pdf)
  #summ
}

# Generate temperature from 25 to 32 degrees, with heavier weight towards the middle
set.seed(123)
Climate_sg <- round(rnorm(366, mean = 28.5, sd = 1.5), 1)
Climate_sg <- pmax(pmin(Climate_sg, 32), 25)  # Ensure within bounds
Temp <- data.frame(Day = seq_along(Climate_sg), Temperature = Climate_sg)
#Temp <- data.frame(Day =  seq_along(Climate_df$temp), Temperature = Climate_df$temp)
Tmax <- length(Temp$Day) # Total number of days
GTmax <- 35              # Maximum generation time (days)

# Define the function for evaluating generation time distribution
evalGenTimeDist <- function(x, a, b, serT, tt, GTmax) {
  mxx <- int_sum_gamma_T(1, a, b, 
                         Temp = serT[x:(Tmax + GTmax + 1)],
                         t = tt[x:(Tmax + GTmax + 1)],       # Time series in days
                         max = GTmax  # Max generation time in days
  )
  return(mxx$pdf)  # Return only the PDF (as in the original code)
}

# Apply the function over the entire time series using mcmapply for parallelization
gt <- mcmapply(
  evalGenTimeDist, 
  1:Tmax, 
  MoreArgs = list(
    a = c(16, 4.3, 1, 1), 
    b = c(1 / 2.69821, 1 / 0.4623722, 1, 1),
    serT = Temp$Temperature,   # Your temperature time series (daily data)
    tt = Temp$Day,         # Your daily time series
    GTmax = GTmax     # Maximum generation time in days
  ), 
  mc.cores = detectCores()  # Use parallel computation with available cores
)

# Load generation time distributions after running previous code only once
load('./Generation_times.RData')

all_weather_data <- read.csv("~/Documents/NTU/Projects/sg_temp_2019.csv")
Temperature <- all_weather_data$Mean.Temperature...C.

num_simulations <- 100
total_cases <- 10000
sim_days <- 366
n_cores <- max(10, (detectCores()) - 1, na.rm = TRUE)
RNGkind("L'Ecuyer-CMRG") # reproducible random numbers with mclapply
set.seed(2026)

# Keep only points that are within the Singapore land area
singapore_sf <- ne_countries(scale = "large", returnclass = "sf") %>%
  filter(sovereignt == "Singapore")

# Function with time-dependent generation interval in sim and calculation ------------------
sim_func_td_td <- function(sim_index, output_directory) {
  # Generate random coordinates with more weight in the center
  lon_min <- 103.70
  lon_max <- 104.03
  lat_min <- 1.23
  lat_max <- 1.47
  
  # Midpoints for the bounds
  lon_center <- (lon_min + lon_max) / 2
  lat_center <- (lat_min + lat_max) / 2
  
  # Standard deviations for concentration around the center
  lon_sd <- (lon_max - lon_min) / 6  # Approx. 99.7% within bounds
  lat_sd <- (lat_max - lat_min) / 6
  
  # Function to generate points with truncated normal distribution
  generate_weighted_points <- function(n) {
    lon <- rtruncnorm(n, a = lon_min, b = lon_max, mean = lon_center, sd = lon_sd)
    lat <- rtruncnorm(n, a = lat_min, b = lat_max, mean = lat_center, sd = lat_sd)
    points <- data.frame(Longitude = lon, Latitude = lat)
    st_as_sf(points, coords = c("Longitude", "Latitude"), crs = 4326)
  }
  
  # Generate coordinates
  candidate_points <- generate_weighted_points(2 * total_cases)
  
  # Keep only points that are within the Singapore land area
  valid_points <- candidate_points %>%
    st_filter(singapore_sf, .predicate = st_within)
  
  # Select the exact number of points needed
  final_points <- valid_points %>%
    slice_sample(n = total_cases)
  
  # Extract coordinates into sim
  coords <- cbind(
    st_coordinates(final_points)[, 1],  # Longitude
    st_coordinates(final_points)[, 2]   # Latitude
  )
  
  # Compute the mean latitude for scaling
  lat_mean <- mean(coords[,2])  # Assuming Longitude is column 1 and Latitude is column 2
  
  # Conversion factors from degrees to meters
  meters_per_degree_lat <- 111132  # 1° latitude ≈ 111.132 km
  meters_per_degree_lon <- 111132 * cos(lat_mean * pi / 180)  # Adjust for latitude
  
  # Convert coordinates to meters
  coords_meters <- coords
  coords_meters[,1] <- coords[,1] * meters_per_degree_lon  # Scale longitude
  coords_meters[,2] <- coords[,2] * meters_per_degree_lat  # Scale latitude
  
  # Compute Euclidean distances in meters
  distance_matrix <- as.matrix(dist(coords_meters, method = "euclidean"))
  
  # Initialize infection indicators, times, and R
  infection_indicator <- rep(NA, total_cases)
  infection_time <- rep(NA, total_cases)
  R_obs <- numeric(total_cases)
  R_exp <- numeric(total_cases)
  
  # Initial index case
  index_case <- sample(1:total_cases, 1)
  infection_indicator[index_case] <- 1
  infection_time[index_case] <- 1
  
  # Initialize list to track interactions
  interaction_list <- list()
  
  # Simulation
  for (t in 1:sim_days) {
    active_indices <- which(!is.na(infection_time) & (t - infection_time) <= 36)
    if (length(active_indices) == 0) next
    
    # Parallelize over active indices (l)
    results <- mclapply(active_indices, function(l) {
      probs <- numeric(total_cases) # vector to hold individual transmission likelihoods
      normalization_factor <- 0
      
      for (k in 1:total_cases) {
        if (!is.na(infection_time[k])) next
        
        dist <- distance_matrix[l, k]
        if (dist > 1500) next
        
        time_diff <- t - infection_time[l]
        if (time_diff < 8 || time_diff > 36) next
        
        temp_index <- infection_time[l]
        Temp <- Temperature[temp_index]
        
        time_weight <- gt[time_diff, Temp]
        spatial_weight <- .99 * exp(- (dist / 106)^2)
        
        weight <- spatial_weight * time_weight
        probs[k] <- weight
        normalization_factor <- normalization_factor + weight
      }
      
      if (normalization_factor > 0) {
        probs <- probs / normalization_factor
      }
      
      list(probs = probs, l = l)
    }, mc.cores = max(23, (detectCores()) - 1, na.rm = TRUE))
    
    # Process results for active indices
    for (k in 1:total_cases) {
      if (!is.na(infection_indicator[k])) next # skip cases that are already infected
      
      # Collect probabilities for all active infectors targeting case k
      infectors <- integer(0)
      infector_probs <- numeric(0)
      
      for (i in results) {
        l <- i$l
        probs <- i$probs
        
        if (probs[k] > 0 && rbern(1, probs[k])) {  # Bernoulli trial for transmission
          infectors <- c(infectors, l)
          infector_probs <- c(infector_probs, probs[k])
        }
      }
      
      if (length(infector_probs) > 0) {
        # Normalize probabilities across all successful infectors for case k
        infector_probs <- infector_probs / sum(infector_probs)
        
        # Draw one infector based on the normalized probabilities
        if (length(infectors) > 1) {
          
          # sample one infector using normalized probabilities
          chosen_infector <- infectors[which.max(cumsum(infector_probs) > runif(1))]
        } else {
          # If only one infector, directly assign it
          chosen_infector <- infectors[1]
        }
        
        # Update the infection and time for case k
        infection_indicator[k] <- 1
        infection_time[k] <- t
        
        # Record the infection event for tracking who infected whom
        interaction_list[[length(interaction_list) + 1]] <-
          list(infector = chosen_infector, infectee = k,
               prob = infector_probs[which.max(cumsum(infector_probs) > runif(1))])
        
        # Update the observed and expected R for the chosen infector
        R_obs[chosen_infector] <- R_obs[chosen_infector] + 1
      }
    }
  }
  
  
  # Create simulation data and calculate MAE
  sim_data <- data.frame(coords) %>%
    mutate(
      ID = 1:total_cases,
      InfectionTime = infection_time,
      R_obs = R_obs,
    ) %>%
    dplyr::select(ID, everything())
  
  calculate_R_exp <- function(data, distance_matrix, gt, mean_dist) {
    # Calculate time differences
    time_diff_matrix <- outer(data$InfectionTime, data$InfectionTime, FUN = function(x, y) as.numeric(x - y))
    
    # Compute spatial weights (exponential decay with mean distance)
    space_weight <- ifelse(distance_matrix <= 1500, .99 * exp(-(distance_matrix/106)^2), 0)
    time_diff_mask <- (time_diff_matrix >= 8) & (time_diff_matrix <= 36)
    space_weight <- space_weight * time_diff_mask
    # Initialize space-time weight matrix
    space_time_weight <- matrix(0, nrow = total_cases, ncol = total_cases)
    
    # Loop through each pair of cases to compute space-time weights
    for (i in 1:total_cases) {
      for (j in 1:total_cases) {
        # Time difference between case i and case j
        time_diff <- time_diff_matrix[i, j]
        
        # Skip invalid time differences
        if (is.na(time_diff) || time_diff < 8 || time_diff > 36) next
        
        # Onset day of case i (used to index the generation time matrix)
        onset_day_i <- data$InfectionTime[i]
        
        # Ensure time_diff and onset_day_i are within bounds of the generation time matrix
        if (time_diff <= nrow(gt) && onset_day_i <= ncol(gt)) {
          # Retrieve the corresponding generation time weight
          gt_value <- gt[time_diff, Temperature[onset_day_i]]
          
          # Skip if the generation time value is NA
          if (!is.na(gt_value)) {
            # Compute space-time weight for this pair
            space_time_weight[i, j] <- gt_value * space_weight[i, j]
          }
        }
      }
    }
    
    # Normalize by all other potential pairs, setting rows with zero sum to zero
    p_ij <- apply(space_time_weight, 1, function(x) {
      row_sum <- sum(x, na.rm = TRUE)
      if (row_sum == 0) {
        return(rep(0, length(x)))  # Set the row to zero if the sum is zero
      } else {
        return(x / row_sum)
      }
    })
    
    num_p_per_l <- rowSums(p_ij > 0)
    num_p_summary <- summary(num_p_per_l)
    var_Rl <- rowSums(p_ij * (1 - p_ij), na.rm = TRUE)
    var_summary <- summary(var_Rl)
    nonzero_p <- p_ij[p_ij > 0]
    p_summary <- summary(nonzero_p)
    max_p_per_l <- apply(p_ij, 1, max)
    max_p_summary <- summary(max_p_per_l)
    
    # Calculate R_exp as the row sums of the normalized probabilities
    R_estim <- rowSums(p_ij, na.rm = TRUE)
    R_estim_df <- R_estim
    return(list(
      R_estim = R_estim_df,
      num_p_per_l = num_p_per_l,
      var_Rl = var_Rl,
      p_summary = p_summary,
      num_p_summary = num_p_summary,
      var_summary = var_summary,
      max_p_summary = max_p_summary
    ))
  }
  
  results_obj <- calculate_R_exp(sim_data, distance_matrix, gt, mean_dist)
  sim_data$R_estim <- results_obj$R_estim
  sim_data$var_Rl <- results_obj$var_Rl
  num_p_summary <- results_obj$num_p_summary
  p_summary <- results_obj$p_summary
  max_p_summary <- results_obj$max_p_summary
  var_summary <- results_obj$var_summary
  sim_data <- sim_data %>% mutate(Temperature = Temperature[InfectionTime])
  daily_theoretical_var <- sim_data %>%
    group_by(InfectionTime) %>%
    summarise(
      R_t_est = mean(R_estim, na.rm = TRUE),
      var_Rt_theoretical = sum(var_Rl, na.rm = TRUE) / n()^2,
      n_infectors = n()
    )
  
  # Bootstrapping within the simulation (parallelized)
  n_bootstraps <- 100  # Number of bootstrap samples
  
  # Define a function to process a single bootstrap sample
  process_bootstrap <- function(b) {
    tryCatch({
      # Resample individuals with replacement
      bootstrap_sample <- sim_data[sample(nrow(sim_data), replace = TRUE), ]
      
      # Recalculate R_estim for the bootstrap sample
      bootstrap_results_obj <- calculate_R_exp(
        bootstrap_sample,
        distance_matrix,
        gt,
        mean_dist
      )
      
      bootstrap_R_estim <- bootstrap_results_obj$R_estim
      # Append R_estim to the bootstrap sample
      bootstrap_sample$R_estim <- bootstrap_R_estim
      
      # Calculate daily R_t for the bootstrap sample
      daily_Rt <- bootstrap_sample %>%
        group_by(InfectionTime) %>%
        summarise(R_t = mean(R_estim, na.rm = TRUE)) %>%
        mutate(Bootstrap = b)  # Add bootstrap index for tracking
      
      return(daily_Rt)
    }, error = function(e) {
      # Log the error message and bootstrap index
      message("Error in bootstrap ", b, ": ", e$message)
      return(NULL)  # Return NULL for failed bootstraps
    })
  }
  
  # Parallelize the bootstrapping loop using mclapply
  bootstrap_Rt <- mclapply(1:n_bootstraps, process_bootstrap, mc.cores = n_cores)
  
  # Remove NULL results (failed bootstraps)
  bootstrap_Rt <- bootstrap_Rt[!sapply(bootstrap_Rt, is.null)]
  
  # Combine all bootstrap results into a single df
  bootstrap_Rt <- do.call(rbind, bootstrap_Rt)
  
  # Check if the bootstrap results are valid
  if (is.null(bootstrap_Rt) || nrow(bootstrap_Rt) == 0) {
    stop("Bootstrap results are empty. Check for errors in the process_bootstrap function.")
  }
  
  # Save the bootstrapped R_t results as a df within R
  bootstrap_results <- bootstrap_Rt
  
  mae <- mean(abs(sim_data$R_obs - round(sim_data$R_estim)), na.rm = TRUE)
  
  bootstrap_output_file <- file.path(output_dir, paste0("bootstrap_results_", sim_index, ".csv"))
  write.csv(bootstrap_Rt, file = bootstrap_output_file, row.names = FALSE)
  output_file <- file.path(output_dir, paste0("simulation_", sim_index, ".csv"))
  write.csv(sim_data, file = output_file, row.names = FALSE)
  
  # Remove large objects explicitly before function ends
  rm(distance_matrix, infection_indicator, infection_time, R_obs, coords)
  
  return(list(
    mae = mae,
    bootstrap_results = bootstrap_results,
    num_p_summary = num_p_summary,
    p_summary = p_summary,
    max_p_summary = max_p_summary,
    var_summary = var_summary,
    daily_theoretical_var = daily_theoretical_var
  ))
}

output_dir <- "."
sim_results_td_td <- lapply(1:num_simulations,
                            function(sim_idx) sim_func_td_td(sim_idx, output_dir))

# Load all bootstrap results from the saved files ---------
bootstrap_files <- list.files(output_dir, pattern = "bootstrap_results_\\d+\\.csv", full.names = TRUE)
bootstrap_tdtd <- bootstrap_files[str_extract(bootstrap_files, "\\d+") %>% as.numeric() %in% 1:100]
bootstrap_tdtd_all <- map_dfr(bootstrap_tdtd, function(file) {
  df <- read.csv(file)
  df$simulation_id <- str_extract(file, "\\d+") %>% as.numeric()
  df
})

# Get confidence intervals
ci_tdtd <- bootstrap_tdtd_all %>% group_by(InfectionTime) %>%
  summarize(R_t_mean = mean(R_t, na.rm = TRUE),
            R_t_lower = quantile(R_t, 0.025, na.rm = TRUE),
            R_t_upper = quantile(R_t, 0.975, na.rm = TRUE))

# Sensitivity, specificity and AUC ---------
tdtd_classif_Rt <- data.frame(Simulation = integer(), Sensitivity = numeric(),
                              Specificity = numeric(), AUC = numeric(),
                              stringsAsFactors = FALSE)

for (i in 1:100) {
  data <- read.csv(file.path(output_dir, paste0("simulation_", i, ".csv")))
  
  # Rt observed
  Rt_obs <- data %>%
    group_by(InfectionTime) %>%
    summarize(Rt_obs = mean(R_obs, na.rm = TRUE), .groups = "drop")
  
  # CI for estimated Rt
  Rt_ci <- bootstrap_tdtd_all %>%
    filter(Bootstrap == i) %>%
    group_by(InfectionTime) %>%
    summarize(R_t_mean = mean(R_t, na.rm = TRUE),
              R_t_lower = quantile(R_t, 0.025, na.rm = TRUE),
              R_t_upper = quantile(R_t, 0.975, na.rm = TRUE),
              .groups = "drop")
  
  # Combine and filter to only confident cases
  combined_data <- full_join(Rt_obs, Rt_ci, by = "InfectionTime") %>%
    filter(R_t_lower > 1 | R_t_upper < 1) %>%  # exclude ambiguous cases
    mutate(Transmission_true = ifelse(Rt_obs >= 1, 1, 0),
           Transmission_pred = ifelse(R_t_lower > 1, 1, 0))  # CI classification
  
  # Skip if only one class is present
  if (length(unique(combined_data$Transmission_true)) > 1 &&
      length(unique(combined_data$Transmission_pred)) > 1) {
    
    conf_matrix <- confusionMatrix(factor(combined_data$Transmission_pred),
                                   factor(combined_data$Transmission_true),
                                   positive = "1")
    
    roc_obj <- roc(combined_data$Transmission_true, combined_data$R_t_mean)
    
    tdtd_classif_Rt <- rbind(tdtd_classif_Rt, data.frame(
      Simulation = i,
      Sensitivity = conf_matrix$byClass["Sensitivity"],
      Specificity = conf_matrix$byClass["Specificity"],
      AUC = auc(roc_obj)))
  }
}

mean(tdtd_classif_Rt$AUC)
range(tdtd_classif_Rt$AUC)
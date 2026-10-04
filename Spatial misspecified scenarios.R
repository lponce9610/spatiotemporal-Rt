library(dplyr)
library(parallel)

load('./Generation_times.RData')
all_weather_data <- read.csv("~/Documents/NTU/HCP_files/sg_temp_2020.csv")
Temperature <- all_weather_data$Mean.Temperature...C.
RNGkind("L'Ecuyer-CMRG") # reproducible random numbers with mclapply
set.seed(2026)

# function combining all ----------
shape1 <- 20^2 / 9^2
scale1 <- 9^2 / 20
shape2 <- 18.2^2 / 6.1^2
scale2 <- 6.1^2 / 18.2
gt_misspecified1 <- list(shape = shape1, scale = scale1)  # First misspecified gt parameters
gt_misspecified2 <- list(shape = shape2, scale = scale2)  # Second misspecified gt parameters

calculate_all_scenarios <- function(data,
                                    coords = data[, c("X1", "X2")],
                                    gt_temp_dependent,  # Temperature-dependent GT matrix
                                    gt_misspecified1,   # First misspecified GT parameters (shape, scale)
                                    gt_misspecified2,   # Second misspecified GT parameters (shape, scale)
                                    Temperature,
                                    max_distances = c(500, 1000, 1500),
                                    max_time_diffs = c(29, 36, 43),
                                    n_bootstraps = 100,
                                    n_cores = parallel::detectCores() - 1) {
  
  # Add case_id if not present
  if (!"case_id" %in% colnames(data)) {
    data$case_id <- 1:nrow(data)
  }
  
  # Precompute distance matrix once (using all original cases)
  lat_mean <- mean(coords[, 2])
  meters_per_degree_lat <- 111132
  meters_per_degree_lon <- 111132 * cos(lat_mean * pi / 180)
  
  coords_meters <- coords
  coords_meters[, 1] <- coords[, 1] * meters_per_degree_lon
  coords_meters[, 2] <- coords[, 2] * meters_per_degree_lat
  distance_matrix <- as.matrix(dist(coords_meters, method = "euclidean"))
  
  # Initialize result lists
  gauss_results <- list()
  
  # Core calculation function
  calculate_R_core <- function(input_data, gt_type, space_kernel, max_dist, max_time) {
    total_cases <- nrow(input_data)
    time_diff_matrix <- outer(input_data$InfectionTime, input_data$InfectionTime,
                              FUN = function(x, y) as.numeric(x - y))
    
    # Handle spatial weights
    if (space_kernel != "NoSpace") {
      # Get subset of distance matrix for current cases
      case_indices <- match(input_data$case_id, data$case_id)
      dist_mat <- distance_matrix[case_indices, case_indices]
      time_diff_mask <- (time_diff_matrix >= 1) & (time_diff_matrix <= max_time)
      
      if (space_kernel == "gauss") {
        space_weight <- ifelse(dist_mat <= max_dist, .99 * exp(-(dist_mat/106)^2), 0)
      }
    }
    
    # Initialize space-time weight matrix
    space_weight <- space_weight * time_diff_mask
    space_time_weight <- matrix(0, nrow = total_cases, ncol = total_cases)
    
    # Compute generation time weights
    if (gt_type == "temp_dependent") {
      for (i in 1:total_cases) {
        for (j in 1:total_cases) {
          time_diff <- time_diff_matrix[i, j]
          if (is.na(time_diff) || time_diff < 1 || time_diff > max_time) next
          
          onset_day_i <- input_data$InfectionTime[i]
          if (time_diff <= nrow(gt_temp_dependent) && onset_day_i <= ncol(gt_temp_dependent)) {
            gt_value <- gt_temp_dependent[time_diff, Temperature[onset_day_i]]
            space_time_weight[i,j] <- gt_value * space_weight[i,j]
          }
        }
      }
    } else {
      # For misspecified GTs
      gt_params <- if (gt_type == "misspecified1") gt_misspecified1 else gt_misspecified2
      gt_weights <- dgamma(time_diff_matrix, shape = gt_params$shape, scale = gt_params$scale)
      gt_weights[time_diff_matrix < 1 | time_diff_matrix > max_time] <- 0
      space_time_weight <- gt_weights * space_weight
    }
    
    # Normalize and calculate R_i
    p_ij <- apply(space_time_weight, 1, function(x) {
      row_sum <- sum(x, na.rm = TRUE)
      if (row_sum == 0) rep(0, length(x)) else x / row_sum
    })
    
    rowSums(p_ij, na.rm = TRUE)
  }
  
  # Bootstrap function
  run_bootstrap <- function(gt_type, space_kernel, max_dist, max_time) {
    process_bootstrap <- function(b) {
      tryCatch({
        sample_b <- data[sample(nrow(data), replace = TRUE), ]
        R_i_b <- calculate_R_core(sample_b, gt_type, space_kernel, max_dist, max_time)
        
        data.frame(
          InfectionTime = sample_b$InfectionTime,
          R_i = R_i_b,
          Bootstrap = b
        ) %>%
          group_by(InfectionTime) %>%
          summarise(R_t = mean(R_i, na.rm = TRUE)) %>%
          mutate(Bootstrap = b)
      }, error = function(e) NULL)
    }
    
    bootstrap_results <- parallel::mclapply(1:n_bootstraps, process_bootstrap, mc.cores = n_cores)
    bootstrap_results <- bootstrap_results[!sapply(bootstrap_results, is.null)]
    
    if (length(bootstrap_results) == 0) return(NULL)
    
    bootstrap_df <- do.call(rbind, bootstrap_results)
    
    Rt_summary <- bootstrap_df %>%
      group_by(InfectionTime) %>%
      summarise(
        Rt_mean = mean(R_t, na.rm = TRUE),
        Rt_lower = quantile(R_t, 0.025, na.rm = TRUE),
        Rt_upper = quantile(R_t, 0.975, na.rm = TRUE)
      )
    
    # Calculate original Rt estimates
    R_i_original <- calculate_R_core(data, gt_type, space_kernel, max_dist, max_time)
    original_Rt <- data.frame(InfectionTime = data$InfectionTime, R_i = R_i_original) %>%
      group_by(InfectionTime) %>%
      summarise(Rt_mean = mean(R_i, na.rm = TRUE))
    
    list(original = original_Rt, bootstrap = Rt_summary)
  }
  
  # Function to generate column names
  make_colname <- function(prefix, max_dist, max_time, gt_type) {
    paste0(prefix, "_", max_dist, "_", max_time, "_", gt_type)
  }
  
  # Define GT types to process
  gt_types <- c("temp_dependent", "misspecified1", "misspecified2")
  gt_names <- c("gt1", "gt2", "gt3")
  
  # Run all scenarios
  for (i in seq_along(gt_types)) {
    gt_type <- gt_types[i]
    gt_name <- gt_names[i]
    
    # Gaussian kernel scenarios
    for (max_dist in max_distances) {
      for (max_time in max_time_diffs) {
        scenario_name <- make_colname("Rt_gauss", max_dist, max_time, gt_name)
        cat("Running scenario:", scenario_name, "\n")
        
        results <- run_bootstrap(gt_type, "gauss", max_dist, max_time)
        if (!is.null(results)) {
          gauss_results[[scenario_name]] <- results$bootstrap %>%
            rename_with(~ paste0(scenario_name, "_", .), .cols = -InfectionTime)
        }
      }
    }
  }
  
  # Combine all results into data frames
  combine_results <- function(result_list) {
    if (length(result_list) == 0) return(data.frame(InfectionTime = unique(data$InfectionTime)))
    
    base_df <- data.frame(
      InfectionTime = unique(data$InfectionTime),
      Temperature = Temperature[unique(data$InfectionTime)]
    )
    
    for (scenario in names(result_list)) {
      base_df <- left_join(base_df, result_list[[scenario]], by = "InfectionTime")
    }
    
    base_df
  }
  
  rm(distance_matrix, coords_meters, calculate_R_core, run_bootstrap)
  gc()
  
  gaussian = combine_results(gauss_results)
}

# run all 100 data sets
results_gaussian <- lapply(1:100, function(i) {
  message("Running data set ", i)
  
  sim_data <- read.csv(paste0("./simulation_", i, ".csv"))
  
  calculate_all_scenarios(
    data = sim_data,
    gt_temp_dependent = gt,
    gt_misspecified1 = gt_misspecified1,
    gt_misspecified2 = gt_misspecified2,
    Temperature = Temperature,
    max_distances = c(1500, 1000, 500),
    max_time_diffs = c(29, 36, 43),
    n_bootstraps = 100,
    n_cores = 12
  )
})
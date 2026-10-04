library(dplyr)
library(parallel)
library(simcausal) # for rbern
library(stringr)
library(caret) # for sensitivity/specificity
library(pROC)

load('./Generation_times.RData')
all_weather_data <- read.csv("./Temperature.csv")
Temperature <- all_weather_data$Mean.Temperature...C.
RNGkind("L'Ecuyer-CMRG") # reproducible random numbers with mclapply
set.seed(2026)

# non spatial sims ------------------------------
sim_func_nospace <- function(sim_index, output_dir) {
  # Initialize infection indicators, times, and R
  infection_indicator <- rep(NA, total_cases)
  infection_time <- rep(NA, total_cases)
  R_obs <- numeric(total_cases)
  
  # Initial index case
  index_case <- sample(1:total_cases, 1)
  infection_indicator[index_case] <- 1
  infection_time[index_case] <- 1
  
  # Track who infected whom
  interaction_list <- list()
  
  # Simulate over days
  for (t in 1:sim_days) {
    active_indices <- which(!is.na(infection_time) & (t - infection_time) <= 36)
    if (length(active_indices) == 0) next
    
    results <- mclapply(active_indices, function(l) {
      probs <- numeric(total_cases)
      normalization_factor <- 0
      
      for (k in 1:total_cases) {
        if (!is.na(infection_time[k])) next
        
        time_diff <- t - infection_time[l]
        if (time_diff < 1 || time_diff > 36) next
        
        temp_index <- infection_time[l]
        Temp <- Temperature[temp_index]
        
        time_weight <- gt[time_diff, Temp]
        weight <- time_weight  # spatial weight removed
        probs[k] <- weight
        normalization_factor <- normalization_factor + weight
      }
      
      if (normalization_factor > 0) {
        probs <- probs / normalization_factor
      }
      
      list(probs = probs, l = l)
    }, mc.cores = max(23, detectCores() - 1))
    
    # Process new infections
    for (k in 1:total_cases) {
      if (!is.na(infection_indicator[k])) next
      
      infectors <- integer(0)
      infector_probs <- numeric(0)
      
      for (i in results) {
        l <- i$l
        probs <- i$probs
        
        if (probs[k] > 0 && rbern(1, probs[k])) {
          infectors <- c(infectors, l)
          infector_probs <- c(infector_probs, probs[k])
        }
      }
      
      if (length(infector_probs) > 0) {
        infector_probs <- infector_probs / sum(infector_probs)
        
        # Draw one infector based on the normalized probabilities
        if (length(infectors) > 1) {
          
          # sample one infector using normalized probabilities
          chosen_infector <- infectors[which.max(cumsum(infector_probs) > runif(1))]
        } else {
          # If only one infector, directly assign it
          chosen_infector <- infectors[1]
        }
        
        infection_indicator[k] <- 1
        infection_time[k] <- t
        
        interaction_list[[length(interaction_list) + 1]] <- list(
          infector = chosen_infector, infectee = k,
          prob = infector_probs[which.max(cumsum(infector_probs) > runif(1))])
        R_obs[chosen_infector] <- R_obs[chosen_infector] + 1
      }
    }
  }
  
  # Create output data
  sim_data <- data.frame(ID = 1:total_cases, InfectionTime = infection_time, R_obs = R_obs)
  sim_data$Temperature <- Temperature[sim_data$InfectionTime]
  
  # Recalculate R_exp (no spatial component)
  calculate_R_exp_nospace <- function(data, gt) {
    time_diff_matrix <- outer(data$InfectionTime, data$InfectionTime, function(x, y) as.numeric(x - y))
    time_diff_mask <- (time_diff_matrix >= 1) & (time_diff_matrix <= 36)
    space_time_weight <- matrix(0, nrow = total_cases, ncol = total_cases)
    
    for (i in 1:total_cases) {
      for (j in 1:total_cases) {
        time_diff <- time_diff_matrix[i, j]
        
        if (is.na(time_diff) || time_diff < 1 || time_diff > 36) next
        
        onset_day_i <- data$InfectionTime[i]
        
        if (onset_day_i > ncol(gt)) next
        
        gt_value <- gt[time_diff, Temperature[onset_day_i]]
        
        if (!is.na(gt_value)) {
          space_time_weight[i, j] <- gt_value
        }
      }
    }
    
    p_ij <- apply(space_time_weight, 1, function(x) {
      row_sum <- sum(x, na.rm = TRUE)
      if (row_sum == 0) return(rep(0, length(x))) else return(x / row_sum)
    })
    
    R_estim <- rowSums(p_ij, na.rm = TRUE)
    
    # Scale rows so none exceed 30
    scaling_factors <- pmin(1, 30 / R_estim)
    p_ij <- p_ij * scaling_factors  # shrink probabilities row-wise
    
    R_estim <- rowSums(p_ij, na.rm = TRUE)
    R_estim_df <- R_estim
    return(R_estim_df)
  }
  
  R_estim_results <- calculate_R_exp_nospace(sim_data, gt)
  sim_data <- cbind(sim_data, R_estim_results)
  sim_data <- sim_data %>% mutate(Temperature = Temperature[InfectionTime])
  
  # Bootstrap R_t estimates
  n_bootstraps <- 100
  process_bootstrap <- function(b) {
    tryCatch({
      bootstrap_sample <- sim_data[sample(nrow(sim_data), replace = TRUE), ]
      bootstrap_sample$R_estim <- calculate_R_exp_nospace(bootstrap_sample, gt)
      bootstrap_sample %>%
        group_by(InfectionTime) %>%
        summarize(R_t = mean(R_estim, na.rm = TRUE)) %>%
        mutate(Bootstrap = b)
    }, error = function(e) {
      message("Error in bootstrap ", b, ": ", e$message)
      return(NULL)
    })
  }
  
  bootstrap_Rt <- mclapply(1:n_bootstraps, process_bootstrap, mc.cores = n_cores)
  bootstrap_Rt <- bootstrap_Rt[!sapply(bootstrap_Rt, is.null)]
  bootstrap_results <- do.call(rbind, bootstrap_Rt)
  
  mae <- mean(abs(sim_data$R_obs - round(sim_data$R_estim)), na.rm = TRUE)
  
  # Save outputs
  write.csv(sim_data, file = file.path(output_dir, paste0("simulation_", sim_index, ".csv")), row.names = FALSE)
  write.csv(bootstrap_results, file = file.path(output_dir, paste0("bootstrap_results_", sim_index, ".csv")), row.names = FALSE)
  
  return(list(mae = mae, bootstrap_results = bootstrap_results))
}

output_dir <- "~/Documents/NTU/Projects/Individual_Rt/No space sims"
num_simulations <- 100
total_cases <- 10000
sim_days <- 366
n_cores <- detectCores() -1
nospace_td_td <- lapply(1:num_simulations,
                        function(sim_index) sim_func_nospace(sim_index, output_dir))


# Amplitude sensitivity results -----------------------
# A = .8 (sims 1-100), A = .5 (sims 101-200), A = .3 (sims 201-300)
bootstrap_files <- list.files(".", pattern = "bootstrap_results_\\d+\\.csv", full.names = TRUE)

amplitude_results <- function(offset) {
  percent_errors <- numeric(100)
  classif_Rt <- data.frame(Simulation = integer(), Sensitivity = numeric(),
                           Specificity = numeric(), AUC = numeric(),
                           stringsAsFactors = FALSE)
  
  # Load all bootstrap results for this amplitude
  bootstrap_set <- bootstrap_files[str_extract(bootstrap_files, "\\d+") %>% as.numeric() %in% (offset + 1:100)]
  bootstrap_all <- lapply(bootstrap_set, read.csv) %>% bind_rows()
  
  for (i in 1:100) {
    data <- read.csv(paste0("./simulation_", i + offset, ".csv"))
    
    # Compute Rt observed and estimated per InfectionTime
    Rt_obs <- data %>% group_by(InfectionTime) %>% summarize(Rt_obs = mean(R_obs, na.rm = TRUE))
    Rt_estim <- data %>% group_by(InfectionTime) %>%
      summarize(Rt_estim = mean(R_estim_results, na.rm = TRUE))
    Daily_Rt <- full_join(Rt_obs, Rt_estim, by = "InfectionTime")
    
    # Compute MAE and percent error
    mae_Rt <- mean(abs(Daily_Rt$Rt_obs - round(Daily_Rt$Rt_estim)), na.rm = TRUE)
    percent_errors[i] <- round((mae_Rt / mean(Daily_Rt$Rt_obs, na.rm = TRUE)) * 100, 2)
    
    # CI for estimated Rt
    Rt_ci <- bootstrap_all %>%
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
      
      classif_Rt <- rbind(classif_Rt, data.frame(
        Simulation = i,
        Sensitivity = conf_matrix$byClass["Sensitivity"],
        Specificity = conf_matrix$byClass["Specificity"],
        AUC = auc(roc_obj)))
    }
  }
  
  data.frame(percent_errors) %>%
    summarize(PercentError_mean = mean(percent_errors, na.rm = TRUE),
              PercentError_min = min(percent_errors, na.rm = TRUE),
              PercentError_Q1 = quantile(percent_errors, 0.25, na.rm = TRUE),
              PercentError_Q3 = quantile(percent_errors, 0.75, na.rm = TRUE),
              PercentError_max = max(percent_errors, na.rm = TRUE)) %>%
    cbind(classif_Rt %>%
            summarize(AUC_mean = mean(AUC, na.rm = TRUE),
                      AUC_min = min(AUC, na.rm = TRUE),
                      AUC_max = max(AUC, na.rm = TRUE)))
}

A80_errors_auc <- amplitude_results(0)
A50_errors_auc <- amplitude_results(100)
A30_errors_auc <- amplitude_results(200)
amplitude_errors_auc <- A80_errors_auc %>% rbind(A50_errors_auc) %>% rbind(A30_errors_auc)
library(tidyverse)
library(fpp2)
library(tseries)
library(urca)
library(xts)
library(ggfortify)
library(Metrics)

# 1. DATA LOADING & PREPARATION

# Data loading
data_df <- read.csv("ts_electric_price.csv", sep = " ")

ts_prices <- ts(data_df[, c(3, 4)], frequency = 365.25)
colnames(ts_prices) <- c("minimum", "maximum")


# 2. UNIVERSE OF DISCOURSE & FUZZY SETS SETUP

# Universe of discourse with a 1% safety margin
D_global <- (max(ts_prices[, 2]) - min(ts_prices[, 1])) * 0.01
Univ <- c(min(ts_prices[, 1]) - D_global, max(ts_prices[, 2]) + D_global)

# Definition of the 20 interval fuzzy sets A1,...,A20
A <- list(
  c(Univ[1], 40), c(0, 60), c(5, 65), c(10, 70), c(20, 80),
  c(30, 90), c(40, 100), c(45, 105), c(50, 110), c(55, 115),
  c(60, 120), c(65, 125), c(70, 130), c(80, 140), c(90, 150),
  c(110, 170), c(130, 190), c(150, 210), c(170, 230), c(190, Univ[2])
)

# 3. SIMILARITY & VALIDATION PARAMETERS

source("SIMILARITIES.R") # Loads the Similarity(x, y, E, Tn, U) function

E_val <- 2  
Tn_val <- 1 
U_val <- 1   

horizon <- 3         # 3-day forecasting horizon
eval_start <- 1067   # Cross-validation start index
eval_end <- nrow(ts_prices) - horizon
steps <- seq(eval_start, eval_end, by = 3) # Step size matches horizon


# 4. CROSS-VALIDATION LOOP (Rolling Forecast)

final_results <- data.frame()
cat("Starting Cross-Validation: 5 FTS Methods vs ARIMA...\n")

for (v in steps) {
  cat("Evaluating fold with training up to day:", v, "\n")
  
  train_v <- ts_prices[1:v, c(1, 2)] 
  actual_v <- ts_prices[(v + 1):(v + horizon), c(1, 2)]
  
  # --- A. ARIMA MODELS ---
  fit_ari_min <- auto.arima(train_v[, 1])
  fit_ari_max <- auto.arima(train_v[, 2])
  pred_ari <- cbind(as.numeric(forecast(fit_ari_min, h = horizon)$mean),
                    as.numeric(forecast(fit_ari_max, h = horizon)$mean))
  
  # --- B. FTS TRAINING ---
  train_norm <- (train_v - Univ[1]) / (Univ[2] - Univ[1])
  A_norm <- lapply(A, function(x) (x - Univ[1]) / (Univ[2] - Univ[1]))
  
  # Historical Fuzzification
  FLR <- sapply(1:v, function(i) {
    which.max(sapply(1:length(A), function(j) 
      Similarity(A_norm[[j]], as.numeric(train_norm[i, ]), E_val, Tn_val, U_val)))
  })
  
  # Fuzzy Logical Relationship Groups (GFLR) & Frequencies
  GFLR <- vector("list", length(A))
  for (i in 1:(v - 1)) GFLR[[FLR[i]]] <- c(GFLR[[FLR[i]]], FLR[i+1])
  GFLR_unique <- lapply(GFLR, unique)
  
  frequencies <- lapply(GFLR, function(tr) {
    if (length(tr) == 0) return(rep(0, length(A)))
    as.numeric(prop.table(table(factor(tr, levels = 1:length(A)))))
  })
  
  # --- C. RECURSIVE FORECASTING (3 DAYS) ---
  p_chen <- list(); p_cheng <- list(); p_sim1 <- list(); p_sim2 <- list(); p_sim3 <- list()
  state_chen <- state_cheng <- state_sim1 <- state_sim2 <- state_sim3 <- FLR[v]
  last_v_chen <- last_v_cheng <- last_v_sim1 <- last_v_sim2 <- last_v_sim3 <- as.numeric(train_v[v, ])
  
  for (h in 1:horizon) {
    
    # Internal function to calculate weights and predict
    predict_v <- function(state, method_type, last_val) {
      next_states <- sort(as.numeric(GFLR_unique[[state]]))
      
      # Persistence rule for unmapped states
      if (length(next_states) == 0) return(as.numeric(train_v[v, ])) 
      
      e_inf <- sapply(next_states, function(s) A[[s]][1])
      e_sup <- sapply(next_states, function(s) A[[s]][2])
      l_norm <- (last_val - Univ[1]) / (Univ[2] - Univ[1])
      
      if (method_type == "chen")  { weights <- rep(1/length(next_states), length(next_states)) }
      if (method_type == "cheng") { weights <- frequencies[[state]][next_states] }
      if (method_type == "sim1")  { 
        w <- sapply(next_states, function(s) Similarity(A_norm[[s]], l_norm, E_val, Tn_val, U_val))
        weights <- if(sum(w) == 0) rep(1/length(next_states), length(next_states)) else w / sum(w)
      }
      if (method_type == "sim2")  { 
        w_s <- sapply(next_states, function(s) Similarity(A_norm[[s]], l_norm, E_val, Tn_val, U_val))
        w_s <- if(sum(w_s) == 0) rep(1/length(next_states), length(next_states)) else w_s / sum(w_s)
        weights <- (w_s + frequencies[[state]][next_states]) / 2
      }
      if (method_type == "sim3")  { 
        w_s <- sapply(next_states, function(s) Similarity(A_norm[[s]], l_norm, E_val, Tn_val, U_val))
        w_s <- if(sum(w_s) == 0) rep(1/length(next_states), length(next_states)) else w_s / sum(w_s)
        weights <- exp((log(w_s + 1e-10) + log(frequencies[[state]][next_states] + 1e-10)) / 2)
      }
      weights <- weights / sum(weights)
      return(c(sum(e_inf * weights), sum(e_sup * weights)))
    }
    
    # Run predictions for each method
    p_chen[[h]]  <- predict_v(state_chen, "chen", last_v_chen)
    p_cheng[[h]] <- predict_v(state_cheng, "cheng", last_v_cheng)
    p_sim1[[h]]  <- predict_v(state_sim1, "sim1", last_v_sim1)
    p_sim2[[h]]  <- predict_v(state_sim2, "sim2", last_v_sim2)
    p_sim3[[h]]  <- predict_v(state_sim3, "sim3", last_v_sim3)
    
    # Update state history for the next step in the horizon
    update_state <- function(val) which.max(sapply(1:length(A), function(j) 
      Similarity(A_norm[[j]], (val - Univ[1]) / (Univ[2] - Univ[1]), E_val, Tn_val, U_val)))
    
    state_chen  <- update_state(p_chen[[h]]);  last_v_chen  <- p_chen[[h]]
    state_cheng <- update_state(p_cheng[[h]]); last_v_cheng <- p_cheng[[h]]
    state_sim1  <- update_state(p_sim1[[h]]);  last_v_sim1  <- p_sim1[[h]]
    state_sim2  <- update_state(p_sim2[[h]]);  last_v_sim2  <- p_sim2[[h]]
    state_sim3  <- update_state(p_sim3[[h]]);  last_v_sim3  <- p_sim3[[h]]
  }
  
  # --- D. FOLD METRIC CALCULATION (Average Similarity) ---
  calc_sim_fold <- function(plist) {
    p_matrix <- do.call(rbind, plist)
    mean(sapply(1:horizon, function(i) {
      Similarity((as.numeric(actual_v[i, ]) - Univ[1]) / (Univ[2] - Univ[1]), 
                 (as.numeric(p_matrix[i, ]) - Univ[1]) / (Univ[2] - Univ[1]), 
                 1, 2, 1) # Metric evaluation using Length and Minimum
    }))
  }
  
  final_results <- rbind(final_results, data.frame(
    Origin = v,
    Chen = calc_sim_fold(p_chen),
    Cheng = calc_sim_fold(p_cheng),
    Proposed_Sim1 = calc_sim_fold(p_sim1),
    Proposed_Sim2 = calc_sim_fold(p_sim2),
    Proposed_Sim3 = calc_sim_fold(p_sim3),
    ARIMA = mean(sapply(1:horizon, function(i) {
      Similarity((as.numeric(actual_v[i, ]) - Univ[1]) / (Univ[2] - Univ[1]), 
                 (as.numeric(pred_ari[i, ]) - Univ[1]) / (Univ[2] - Univ[1]), 
                 1, 2, 1)
    }))
  ))
}


# 6. FINAL RESULTS & VISUALIZATION

cat("\n--- FINAL COMPARISON (AVERAGE SIMILARITY) ---\n")
final_table <- colMeans(final_results[, -1])
print(round(final_table, 4))

# Plot similarity evolution across cross-validation origins
results_long <- final_results %>% pivot_longer(-Origin, names_to = "Method", values_to = "Similarity")

ggplot(results_long, aes(x = Origin, y = Similarity, color = Method)) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  labs(
    title = "Average similarity evolution by method", 
    subtitle = paste("Horizon:", horizon, "days | Training Similarity Measure E =", E_val),
    x = "Forecast Origin Day", 
    y = "Similarity Index"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


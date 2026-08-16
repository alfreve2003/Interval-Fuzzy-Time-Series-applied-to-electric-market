library(tidyverse)
library(fpp2)
library(tseries)
library(urca)
library(xts)
library(ggfortify)
library(Metrics)

# 1. DATA LOADING & PREPARATION

data_df <- read.csv("ts_electric_price.csv", sep = " ")
data_df$date <- as.Date(data_df$fecha)

# Create time series using Minimum Price and Maximum Price
ts_prices <- ts(data_df[, c(3, 4)], frequency = 365.25)
colnames(ts_prices) <- c("minimum_price", "maximum_price")

# Train-Test Split (80-20)
train_data <- ts_prices[1:1115, ]
test_data  <- ts_prices[1116:1118, ]

train_len <- nrow(train_data)
test_len  <- nrow(test_data)

# 2. ARIMA MODELS

fit_ari_min <- auto.arima(train_data[, 1])
fit_ari_max <- auto.arima(train_data[, 2])
pred_arima  <- cbind(
  as.numeric(forecast(fit_ari_min, h = test_len)$mean),
  as.numeric(forecast(fit_ari_max, h = test_len)$mean)
)


# 3. FUZZY TIME SERIES

# Universe of Discourse with 1% safety margin
D_global <- (max(train_data[, 2]) - min(train_data[, 1])) * 0.01
Univ <- c(min(train_data[, 1]) - D_global, max(train_data[, 2]) + D_global)

# EDA
summary(data_df$precio_minimo)
summary(data_df$precio_maximo)

# Plot price distribution histograms
par(mfrow = c(1, 2))
hist(data_df$precio_minimo, main = "Minimum Price", col = "skyblue", xlab = "Price (EUR/MWh)", ylab = "Frequency")
hist(data_df$precio_maximo, main = "Maximum Price", col = "salmon", xlab = "Price (EUR/MWh)", ylab = "Frequency")
par(mfrow = c(1, 1))

summary(data_df$precio_maximo - data_df$precio_minimo)
quantile(data_df$precio_minimo)
quantile(data_df$precio_maximo)

# Definition of 20 Fuzzy Sets
A <- list(
  c(Univ[1], 40), c(0, 60), c(5, 65), c(10, 70), c(20, 80),
  c(30, 90), c(40, 100), c(45, 105), c(50, 110), c(55, 115),
  c(60, 120), c(65, 125), c(70, 130), c(80, 140), c(90, 150),
  c(110, 170), c(130, 190), c(150, 210), c(170, 230), c(190, Univ[2])
)

levels_A <- paste0("A", 1:20)

# Normalization
train_norm <- (train_data - Univ[1]) / (Univ[2] - Univ[1])
A_norm     <- lapply(A, function(j) (j - Univ[1]) / (Univ[2] - Univ[1]))


#  4. FTS TRAINING & FUZZIFICATION

source("SIMILARITIES.R") # Loads Similarity(x, y, E, Tn, U) function

E_val <- 1   # Embedding selection 
# 1: Length, 2: Lukasiewicz, 3: Godel, 4: Goguen, 5: Fodor, 6: Exponential
Tn_val <- 1  # Aggregation function 
# 1: Arithmetic mean, 2: Minimum, 3: Geometric mean, 4: Maximum
U_val <- 1   # Aggregated similarity (1: Yes, 0: No)

# Fuzzification 
FLR <- sapply(1:train_len, function(i) {
  which.max(sapply(1:length(A), function(j) {
    Similarity(A_norm[[j]], as.numeric(train_norm[i, ]), E_val, Tn_val, U_val)
  }))
})

# Fuzzy Logical Relationship Groups (FLRG)
GFLR <- vector("list", length(A))
names(GFLR) <- paste0("A", 1:length(A))

for (i in 1:(train_len - 1)) {
  GFLR[[FLR[i]]] <- c(GFLR[[FLR[i]]], FLR[i + 1])
}

GFLR_unique <- lapply(GFLR, unique)

frequencies <- lapply(GFLR, function(tr) {
  if (length(tr) == 0) return(rep(0, length(A)))
  as.numeric(prop.table(table(factor(tr, levels = 1:length(A)))))
})


# 5. FORECASTING (DEFFUZYFICATION)

pred_chen  <- vector("list", test_len)
pred_cheng <- vector("list", test_len)
pred_sim1  <- vector("list", test_len)
pred_sim2  <- vector("list", test_len)
pred_sim3  <- vector("list", test_len)

state_chen <- state_cheng <- state_sim1 <- state_sim2 <- state_sim3 <- FLR[train_len]
last_chen  <- last_cheng  <- last_sim1  <- last_sim2  <- last_sim3  <- as.numeric(train_data[train_len, ])

for (h in 1:test_len) {
  # Internal function to calculate weights and predict
  predict_v <- function(state, method_type, last_val) {
    next_states <- sort(as.numeric(GFLR_unique[[state]]))
    
    # Persistence rule if no historical transition exists
    if (length(next_states) == 0) return(as.numeric(train_data[train_len, ]))
    
    e_inf  <- sapply(next_states, function(s) A[[s]][1])
    e_sup  <- sapply(next_states, function(s) A[[s]][2])
    l_norm <- (last_val - Univ[1]) / (Univ[2] - Univ[1])
    
    if (method_type == "chen") {
      weights <- rep(1 / length(next_states), length(next_states))
    }
    if (method_type == "cheng") {
      weights <- frequencies[[state]][next_states]
    }
    if (method_type == "sim1") {
      w <- sapply(next_states, function(s) Similarity(A_norm[[s]], l_norm, E_val, Tn_val, U_val))
      weights <- if (sum(w) == 0) rep(1 / length(next_states), length(next_states)) else w / sum(w)
    }
    if (method_type == "sim2") {
      w_s <- sapply(next_states, function(s) Similarity(A_norm[[s]], l_norm, E_val, Tn_val, U_val))
      w_s <- if (sum(w_s) == 0) rep(1 / length(next_states), length(next_states)) else w_s / sum(w_s)
      weights <- (w_s + frequencies[[state]][next_states]) / 2
    }
    if (method_type == "sim3") {
      w_s <- sapply(next_states, function(s) Similarity(A_norm[[s]], l_norm, E_val, Tn_val, U_val))
      w_s <- if (sum(w_s) == 0) rep(1 / length(next_states), length(next_states)) else w_s / sum(w_s)
      weights <- exp((log(w_s + 1e-10) + log(frequencies[[state]][next_states] + 1e-10)) / 2)
    }
    
    weights <- weights / sum(weights) # Re-normalize weights
    return(c(sum(e_inf * weights), sum(e_sup * weights)))
  }
  
  # Run predictions for each method
  pred_chen[[h]]  <- predict_v(state_chen, "chen", last_chen)
  pred_cheng[[h]] <- predict_v(state_cheng, "cheng", last_cheng)
  pred_sim1[[h]]  <- predict_v(state_sim1, "sim1", last_sim1)
  pred_sim2[[h]]  <- predict_v(state_sim2, "sim2", last_sim2)
  pred_sim3[[h]]  <- predict_v(state_sim3, "sim3", last_sim3)
  
  # Update state and last value for next horizon step
  update_state <- function(val) {
    which.max(sapply(1:length(A), function(j) {
      Similarity(A_norm[[j]], (val - Univ[1]) / (Univ[2] - Univ[1]), E_val, Tn_val, U_val)
    }))
  }
  
  state_chen  <- update_state(pred_chen[[h]]);  last_chen  <- pred_chen[[h]]
  state_cheng <- update_state(pred_cheng[[h]]); last_cheng <- pred_cheng[[h]]
  state_sim1  <- update_state(pred_sim1[[h]]);  last_sim1  <- pred_sim1[[h]]
  state_sim2  <- update_state(pred_sim2[[h]]);  last_sim2  <- pred_sim2[[h]]
  state_sim3  <- update_state(pred_sim3[[h]]);  last_sim3  <- pred_sim3[[h]]
}


# 6. EVALUATION METRICS (RMSE, MAPE, SIMILARITY)

calc_rmse <- function(plist) {
  p_m <- do.call(rbind, plist)
  return(c(rmse(test_data[, 1], p_m[, 1]), rmse(test_data[, 2], p_m[, 2])))
}

calc_mape <- function(plist) {
  p_m <- do.call(rbind, plist)
  return(c(mape(test_data[, 1], p_m[, 1]), mape(test_data[, 2], p_m[, 2])) * 100)
}

calc_sim_metric <- function(plist) {
  p_m <- do.call(rbind, plist)
  mean(sapply(1:test_len, function(i) {
    Similarity(
      (as.numeric(test_data[i, ]) - Univ[1]) / (Univ[2] - Univ[1]),
      (as.numeric(p_m[i, ]) - Univ[1]) / (Univ[2] - Univ[1]),
      6, 2, 1
    )
  }))
}

# Build comparison tables

# Conversión a lista de vectores por fila (igual que pred_chen, pred_cheng, etc.)
pred_arima <- lapply(1:nrow(pred_arima), function(i) as.numeric(pred_arima[i, ]))

rmse_comparison <- data.frame(
  ARIMA = calc_rmse(pred_arima),
  Chen = calc_rmse(pred_chen),
  Cheng = calc_rmse(pred_cheng),
  Proposed_Sim1 = calc_rmse(pred_sim1),
  Proposed_Sim2 = calc_rmse(pred_sim2),
  Proposed_Sim3 = calc_rmse(pred_sim3)
)
rownames(rmse_comparison) <- c("Lower_Bound", "Upper_Bound")

mape_comparison <- data.frame(
  ARIMA = calc_mape(pred_arima),
  Chen = calc_mape(pred_chen),
  Cheng = calc_mape(pred_cheng),
  Proposed_Sim1 = calc_mape(pred_sim1),
  Proposed_Sim2 = calc_mape(pred_sim2),
  Proposed_Sim3 = calc_mape(pred_sim3)
)
rownames(mape_comparison) <- c("Lower_Bound", "Upper_Bound")

similarity_comparison <- data.frame(
  ARIMA = calc_sim_metric(pred_arima),
  Chen = calc_sim_metric(pred_chen),
  Cheng = calc_sim_metric(pred_cheng),
  Proposed_Sim1 = calc_sim_metric(pred_sim1),
  Proposed_Sim2 = calc_sim_metric(pred_sim2),
  Proposed_Sim3 = calc_sim_metric(pred_sim3)
)

cat("\n--- RMSE COMPARISON ---\n")
print(round(rmse_comparison, 4))

cat("\n--- MAPE (%) COMPARISON ---\n")
print(round(mape_comparison, 4))

cat("\n--- AVERAGE SIMILARITY METRIC ---\n")
print(round(similarity_comparison, 4))


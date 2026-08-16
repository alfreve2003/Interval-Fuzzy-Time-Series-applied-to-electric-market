# Fuzzy Time Series (FTS) vs ARIMA: Electricity Price Interval Forecasting

This repository contains the **practical implementation in R** for my Bachelor's Thesis (*Trabajo Fin de Grado - TFG*): **Embedding-based similarities for interval data**. The theoretical framework of the thesis was specifically developed to construct, formalize, and evaluate these functions, whose application focused on the design of time series forecasting algorithms based on a fuzzy logic approach.

The full thesis document (**written in Spanish**) is included in this repository: `TFG.pdf`.


##  Repository Structure
- `TFG.pdf`: Full Bachelor's Thesis document (in Spanish).
- `ts_electric_price.csv`: Time series dataset with daily min/max electricity prices.
- `Similarities.R`: Core functions for fuzzy similarity metric calculations.
- `ARIMA_and_FTS.R`: Main R script: train/test split, FTS & ARIMA fitting, evaluation.
- `Split_time_series.R`: R script for time series cross-validation.


## Prerequisites & Dependencies

To execute the scripts in R, install the following required packages:
`install.packages(c("tidyverse", "fpp2", "tseries", "urca", "xts", "ggfortify", "Metrics"))`.


## Practical Implementation & Workflow
- Data Preprocessing: loads time series data and defines the Universe of Discourse ($U$) with safety bounds.
- Fuzzification: maps historical observations into 20 interval fuzzy sets ($A_1 \dots A_{20}$) by maximum similarity.
- Fuzzy Logical Relationship Groups (FLRG): establishes transition between fuzzy sets, creating the fuzzy time series.
- Forecasting and Defuzzyfication: computes step-ahead interval predictions across 5 FTS methods (Chen, Cheng, Sim1, Sim2, Sim3, last three were originally proposed in this project).
- Model Evaluation: evaluates precision using RMSE, MAPE, and Average Similarity indices across predicted intervals, and establishes comparisons between these predictions and those obtained with ARIMA.


## Evaluation Metrics
- RMSE & MAPE: calculated independently for upper and lower price bounds.
- Average similarity: similarity measures were used in order to obtain a global-performance evaluation metric.

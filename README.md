# \# Robust Weighted Ridge Regression: Simulation and Real-Data Examples

# 

# This repository provides code to (i) apply \*\*robust weighted ridge regression\*\* to several benchmark datasets and (ii) run a \*\*Monte Carlo simulation\*\* comparing classical and robust ridge-type estimators under \*\*multicollinearity\*\* and \*\*outlier contamination\*\*.

# 

# \## Contents

# \- \*\*Real-life data examples\*\*

# &nbsp; - `bodyfat` (package: `isdals`)

# &nbsp; - `gasoline.xlsx` (Excel file included in this repository)

# &nbsp; - `cement` and `seatpos` (package: `faraway`)

# \- \*\*Simulation study\*\*

# &nbsp; - `ridgeglr.simu(n, p, nn)` Monte Carlo design

# 

# \## Requirements

# R (≥ 4.0 recommended) and the following packages:

# \- `readxl` (for Excel import)

# \- `isdals` (for `bodyfat` dataset)

# \- `faraway` (for `cement`, `seatpos` datasets)

# 

# Install (if needed):

# 

# ```r

# install.packages(c("readxl", "isdals", "faraway"))

# ```

# 

# \## Quick start

# 1\. Place scripts under `src/` and (if needed) source them:

# 

# ```r

# \# Example: source project functions

# \# source("src/robustweightedRidge.R")

# \# source("src/ridgeglr\_simu.R")

# ```

# 

# 2\. Run a real-data example or the simulation (see below).

# 

# \## Real-life data examples

# \### 1) Bodyfat data

# 

# ```r

# install.packages("isdals")

# library(isdals)

# data(bodyfat)

# 

# \# Define X and y

# bodyfat\_x <- bodyfat\[, -1]

# bodyfat\_y <- bodyfat\[,  1]

# 

# bodyfat.ridge <- robustweightedRidge(bodyfat\_x, bodyfat\_y)

# bodyfat.ridge$esttabledataframe

# bodyfat.ridge$stdtabledataframe

# ```

# 

# \### 2) Gasoline data (Excel)

# Download or copy `gasoline.xlsx` into the repository under `data/`.

# 

# ```r

# install.packages("readxl")

# library(readxl)

# 

# \# Read the gasoline dataset (.xlsx) into R

# gasoline <- read\_excel("data/gasoline.xlsx", sheet = 1)

# 

# \# Define X and y

# gasoline\_x <- gasoline\[, 2:5]

# gasoline\_y <- gasoline\[, 1]

# 

# gasoline.ridge <- robustweightedRidge(gasoline\_x, gasoline\_y)

# gasoline.ridge$esttabledataframe

# gasoline.ridge$stdtabledataframe

# ```

# 

# \### 3) Portland Cement data

# 

# ```r

# install.packages("faraway")

# library(faraway)

# data(cement)

# 

# \# Define X and y

# cement\_x <- cement\[, 1:4]

# cement\_y <- cement\[, 5]

# 

# cement.ridge <- robustweightedRidge(cement\_x, cement\_y)

# cement.ridge$esttabledataframe

# cement.ridge$stdtabledataframe

# ```

# 

# \### 4) Seat position data

# 

# ```r

# install.packages("faraway")

# library(faraway)

# data(seatpos)

# 

# \# Define X and y

# seatpos\_x <- as.matrix(seatpos\[, -9])  # 9th column is 'hipcenter' (target)

# seatpos\_y <- seatpos$hipcenter

# 

# seatpos.ridge <- robustweightedRidge(seatpos\_x, seatpos\_y)

# seatpos.ridge$esttabledataframe

# seatpos.ridge$stdtabledataframe

# ```

# 

# \## Simulation study

# \### Purpose

# `ridgeglr.simu` implements a Monte Carlo simulation study to evaluate and compare the performance of ridge regression estimators, with an emphasis on \*\*robustness\*\* and \*\*efficiency\*\* under:

# \- strong \*\*multicollinearity\*\* (correlated predictors), and

# \- \*\*outlier contamination\*\* in the error distribution.

# 

# \### Methods compared

# \- \*\*LS\*\*: Ordinary Least Squares

# \- \*\*LS Ridge\*\*: Ridge Regression (LS-based)

# \- \*\*S-type Ridge\*\*: Robust ridge using S-type estimation

# \- \*\*Huber M-Ridge\*\*: Weighted ridge with Huber’s M-estimator

# \- \*\*Tukey M-Ridge\*\*: Weighted ridge with Tukey’s bisquare M-estimator

# \- \*\*S-Ridge\*\*: Weighted ridge with S-estimator

# \- \*\*MM-Ridge\*\*: Weighted ridge with MM-estimator

# 

# \### Simulation design (high level)

# \- \*\*Multicollinearity\*\*: Predictors are generated using a user-specified correlation structure. The correlation intensity can be modified by choosing alternative correlation/covariance matrices in the script.

# \- \*\*Error distribution / outliers\*\*: Errors follow a contaminated normal model (e.g., 10% contamination with scale factor `k = 3`). For a purely Gaussian setting, the contamination component can be disabled in the script.

# 

# \### Inputs

# \- `n`: sample size

# \- `p`: number of predictors (\*\*fixed at `p = 3`\*\* in this configuration)

# \- `nn`: number of Monte Carlo replications

# 

# \### Outputs

# \- `simutabledataframe`: Mean, Bias, Variance, MSE, and Relative Efficiency (RE) for all estimators across parameters (β₀, β₁, β₂, β₃)

# \- `meanrankdataframe`: average ranks of methods based on MSE and associated standard errors

# 

# \### Example

# 

# ```r

# ridgesimu <- ridgeglr.simu(n = 100, p = 3, nn = 3000)

# ```

# 




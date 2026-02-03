## Real-life data and simulation examples

### Bodyfat data

```r
install.packages("isdals")
library(isdals)
data(bodyfat)
# Define X and y
bodyfat_x <- bodyfat[,-1]
bodyfat_y <- bodyfat[,1]
bodyfat.ridge=robustweightedRidge(bodyfat_x,bodyfat_y)
bodyfat.ridge$esttabledataframe
bodyfat.ridge$stdtabledataframe
```

### Gasoline data

First, the dataset should be downloaded to the local computer and stored in a designated directory. 
Subsequently, the file path of this directory can be specified in R, and the dataset named gasoline can be imported
into the R environment using the appropriate data‐loading command given below.

```r
install.packages("readxl")
library(readxl)

# Read the gasoline dataset (.xlsx) into R
gasoline <- read_excel("gasoline.xlsx", sheet = 1)

# Define X and y
gasoline_x=gasoline[,2:5]
gasoline_y=gasoline[,1]
gasoline.ridge=robustweightedRidge(gasoline_x,gasoline_y)
gasoline.ridge$esttabledataframe
gasoline.ridge$stdtabledataframe
```

### PortlandCement data

```r
install.packages("faraway")
library(faraway)
data(cement)

# Define X and y
cement_x=cement[,1:4]
cement_y=cement[,5]
cement.ridge=robustweightedRidge(cement_x,cement_y)
cement.ridge$esttabledataframe
cement.ridge$stdtabledataframe
```

### seatpos

```r
install.packages("faraway")
library(faraway)
data(seatpos)

# Define X and y
seatpos_x <- as.matrix(seatpos[, -9]) # 9th column is 'hipcenter' (target)
seatpos_y <- seatpos$hipcenter
seatpos.ridge=robustweightedRidge(seatpos_x,seatpos_y)
seatpos.ridge$esttabledataframe
seatpos.ridge$stdtabledataframe
```

### simulation

**Purpose:** This function implements a Monte Carlo simulation study to evaluate and compare the performance of various ridge regression estimators. The simulation is specifically designed to assess estimator robustness and efficiency under conditions of high multicollinearity and outlier contamination in the error distribution.

**Methods Compared:**
- LS: Ordinary Least Squares
- LS Ridge: Ridge Regression with LS
- S-type Ridge: Robust ridge using S-type estimator
- Huber M-Ridge: Weighted ridge with Huber M-estimator
- Tukey M-Ridge: Weighted ridge with Tukey bisquare
- S-Ridge: Weighted ridge with S-estimator
- MM-Ridge: Weighted ridge with MM-estimator

**Simulation Design and Data Generation:**
- **Multicollinearity:** Data are generated with a specified correlation structure among predictors. The correlation intensity can be adjusted by selecting alternative covariance matrices within the script.
- **Error Distribution and Outliers:** To evaluate robustness, the error term is modeled using a contaminated normal distribution (e.g., 10% contamination with a scale factor of k=3). The script allows for a transition to a standard normal distribution by modifying the error model parameters.

**Output Metrics:**
- `simutabledataframe`: Comprehensive performance metrics including Mean, Bias, Variance, Mean Squared Error (MSE), and Relative Efficiency (RE) for all estimators across all parameters (β0, β1, β2, β3)
- `meanrankdataframe`: Statistical ranking of the methods based on MSE performance and associated standard errors

**Input Metrics:**
- `n`  = sample size
- `p`  = Number of predictors (Fixed at p=3 for this configuration)
- `nn` = number of simulation replications

**Usage:**

```r
ridgesimu=ridgeglr.simu(n=100,p=3,nn=3000)
```

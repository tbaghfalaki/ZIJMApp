
```
# Clear all objects from the workspace
rm(list=ls())

# Load the necessary libraries
library(JM)
library(ZI)
library(DPCri)
library(survival)

# Load the AIDS dataset
data("aids")
# Transform CD4 data by squaring and rounding
y = round(aids$CD4^2)
# Display the minimum value of the transformed CD4 data
min(y)
# Add transformed CD4 data to the dataset
aids$Y = y
# Attach the dataset for easier variable access
attach(aids)

# Set seed for reproducibility
set.seed(9)

# Randomly select 80% of the data for the training set
INDTRAIN <- sample(aids.id$patient, round(0.8 * dim(aids.id)[1]))

# Create training longitudinal data
dataLong_t <- subset(aids, aids$patient %in% INDTRAIN)

# Create training survival data
dataSurv_t <- subset(aids.id, aids.id$patient %in% INDTRAIN)

# Display column names of the training survival data
names(dataSurv_t)

# Create validation longitudinal data
dataLong_v <- subset(aids, !(aids$patient %in% INDTRAIN))

# Create validation survival data
dataSurv_v <- subset(aids.id, !(aids.id$patient %in% INDTRAIN))

# Fit a Zero-Inflated Joint Model with specified fixed and random effects
Z1NB <- ZIJMCV(
  FixedY = Y ~ gender + obstime + drug * obstime + prevOI + AZT, 
  RandomY = ~obstime, 
  GroupY = ~patient,
  FixedZ = ~ obstime + drug, 
  RandomZ = ~obstime, 
  GroupZ = ~patient,
  formSurv = survival::Surv(Time, death) ~ gender + drug + prevOI + AZT,
  dataLong = dataLong_t, 
  dataSurv = dataSurv_t,
  obstime = "obstime", 
  id = "patient",
  n.chains = 2,
  n.iter = 200000, 
  n.burnin = 150000, 
  n.thin = 60, 
  K = 15, 
  family = "NB"
)

# Save the fitted model
save(Z1NB, file="Z1_outNB2.RData")

# If needed, load the fitted model
# load("Z1_outNB.RData")

# Initialize parameters for dynamic prediction
s1 = 0
t1 = 3

# Compute dynamic prediction for s1 = 0
DD <- DP_CV(
  object = Z1NB, s = s1, t = t1, 
  n.chains = 1, n.iter = 10000, n.burnin = 5000, 
  n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
)

# Calculate criteria (e.g., AUC, Brier Score) for the prediction
A1 = Criteria(
  s = s1, t = t1, 
  Survt = dataSurv_v$Time, 
  CR = dataSurv_v$death, 
  P = DD$DP$est, 
  cause = 1
)$Cri

# Update s1 and compute dynamic prediction for s1 = 3
s1 = 3
DD <- DP_CV(
  object = Z1NB, s = s1, t = t1, 
  n.chains = 1, n.iter = 10000, n.burnin = 5000, 
  n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
)
A2 = Criteria(
  s = s1, t = t1, 
  Survt = dataSurv_v$Time, 
  CR = dataSurv_v$death, 
  P = DD$DP$est, 
  cause = 1
)$Cri

# Update s1 and compute dynamic prediction for s1 = 6
s1 = 6
DD <- DP_CV(
  object = Z1NB, s = s1, t = t1, 
  n.chains = 1, n.iter = 10000, n.burnin = 5000, 
  n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
)
A3 = Criteria(
  s = s1, t = t1, 
  Survt = dataSurv_v$Time, 
  CR = dataSurv_v$death, 
  P = DD$DP$est, 
  cause = 1
)$Cri

# Update s1 and compute dynamic prediction for s1 = 9
s1 = 9
DD <- DP_CV(
  object = Z1NB, s = s1, t = t1, 
  n.chains = 1, n.iter = 10000, n.burnin = 5000, 
  n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
)
A4 = Criteria(
  s = s1, t = t1, 
  Survt = dataSurv_v$Time, 
  CR = dataSurv_v$death, 
  P = DD$DP$est, 
  cause = 1
)$Cri

# Update s1 and compute dynamic prediction for s1 = 12
s1 = 12
DD <- DP_CV(
  object = Z1NB, s = s1, t = t1, 
  n.chains = 1, n.iter = 10000, n.burnin = 5000, 
  n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
)
A5 = Criteria(
  s = s1, t = t1, 
  Survt = dataSurv_v$Time, 
  CR = dataSurv_v$death, 
  P = DD$DP$est, 
  cause = 1
)$Cri

# Update s1 and compute dynamic prediction for s1 = 15
s1 = 15
DD <- DP_CV(
  object = Z1NB, s = s1, t = t1, 
  n.chains = 1, n.iter = 10000, n.burnin = 5000, 
  n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
)
A6 = Criteria(
  s = s1, t = t1, 
  Survt = dataSurv_v$Time, 
  CR = dataSurv_v$death, 
  P = DD$DP$est, 
  cause = 1
)$Cri
```

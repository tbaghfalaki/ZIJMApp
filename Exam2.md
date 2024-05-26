R code for ZINB joint modeling
---------------

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

```
The output of the package for estimation is as follows:

```
> Z1NB$Estimation
$Count_model
                         Est         SD        L_CI        U_CI     Rhat
(Intercept)      4.314823359 0.22063872  3.85417074  4.73927679 1.071021
gendermale       0.147060241 0.20401965 -0.25589916  0.56513281 1.027733
obstime         -0.082278299 0.01104098 -0.10413032 -0.06090384 1.004544
drugddI         -0.034428255 0.12556266 -0.27998438  0.21260648 1.008524
prevOIAIDS      -1.196801915 0.15427707 -1.46152422 -0.85487660 1.074015
AZTfailure      -0.056102666 0.15663429 -0.35515606  0.25343292 1.094113
obstime:drugddI  0.009411766 0.01565015 -0.01977149  0.04073289 1.000102
Dispersion       6.025612132 0.47599626  5.10965194  6.97086767 1.002437

$Zero_inflated_model
                   Est        SD      L_CI        U_CI     Rhat
(Intercept) -4.2917028 0.9626093 -7.616277 -3.16781112 1.121065
obstime     -0.4383768 0.2939147 -1.120063 -0.05192336 1.009535
drugddI     -0.9990714 0.6548727 -2.391691  0.14064177 1.007103

$Survival_model
                     Est         SD        L_CI       U_CI      Rhat
gendermale   -0.53072775 0.30891769 -1.11249107  0.1268829 1.0003303
drugddI       0.39292943 0.23196269 -0.01255731  0.8768509 1.0006331
prevOIAIDS    1.16057466 0.28195560  0.61569518  1.7130565 1.0058120
AZTfailure    0.26261609 0.21941415 -0.14905035  0.7074987 1.0034592
gamma_lambda -0.38611701 0.07464751 -0.54693312 -0.2464656 0.9996324
gamma_pi      0.10413065 0.10629279 -0.08790489  0.3488038 1.0122079
h1            0.08673818 0.07110144  0.01913615  0.2670304 1.0110370
h2            0.17099632 0.13569690  0.02649546  0.5127538 1.0082621
h3            0.11864850 0.11207866  0.01260471  0.3673736 1.0051442
h4            0.23700319 0.25722658  0.02056355  0.8157685 1.0265779
h5            0.36894298 0.48728188  0.01933185  1.4268546 1.0115499

$D
$D$D11
$D$D11$est
           Intercept      Slope
Intercept 1.32317602 0.01066162
Slope     0.01066162 0.01227506

$D$D11$sd
          Intercept       Slope
Intercept 0.1144397 0.010198001
Slope     0.0101980 0.001381679

$D$D11$L
          Intercept       Slope
Intercept  1.113976 -0.00959100
Slope     -0.009591  0.00974015

$D$D11$U
           Intercept      Slope
Intercept 1.54872870 0.03089683
Slope     0.03089683 0.01520182

$D$D11$Rhat
          Intercept     Slope
Intercept 0.9996436 0.9994549
Slope     0.9994549 1.0002050


$D$D22
$D$D22$est
           Intercept      Slope
Intercept 1.72049097 0.04406276
Slope     0.04406276 0.22672109

$D$D22$sd
          Intercept     Slope
Intercept 2.9608167 0.2768176
Slope     0.2768176 0.2061384

$D$D22$L
           Intercept       Slope
Intercept  0.1640189 -0.37373777
Slope     -0.3737378  0.07042685

$D$D22$U
           Intercept     Slope
Intercept 14.2602526 0.5963163
Slope      0.5963163 0.6714977

$D$D22$Rhat
          Intercept    Slope
Intercept  1.031298 1.085842
Slope      1.085842 1.086394
```

Compute risk prediction, as well as AUC and BS

```
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
The R code for ZIP and ZIGP follows the same structure as described above.

Finally, the individual dynamic prediction (DP) is as follows:
```
DPplot1(Z1NB,
        s = 10, id_new = 184, by = .5, mi = 10,
        Marker_lab="CD4", Time_lab="Time (Month)",
        n.chains = 1, n.iter = 1000, n.burnin = 500,
        dataLong = dataLong_v, dataSurv = dataSurv_v
)
title("# 184, ZINB", font.main=2, cex.main=1)



DPplot1(Z1NB,
        s = 4, id_new = 378, by = .5, mi = 10,
        Marker_lab="CD4", Time_lab="Time (Month)",
        n.chains = 1, n.iter = 1000, n.burnin = 500,
        dataLong = dataLong_v, dataSurv = dataSurv_v
)
title("# 378, ZINB", font.main=2, cex.main=1)

```

resulting in the following plot:


<img src="/Figure/ind_aids.png" alt="Description" width="600" height="500">











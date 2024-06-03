R code for ZINB joint modeling for Pregnency data
---------------

```
# Clear the environment
rm(list=ls())

# Load necessary libraries
library(ZI) # devtools::install_github("tbaghfalaki/ZI")
library(DPCri) # devtools::install_github("tbaghfalaki/DPCri")
library(survival)

# Read in the data
Long = read.table("Long.txt", header=TRUE)
Surv = read.table("surv.txt", header=TRUE)

# Display column names of the datasets
names(Long)
names(Surv)

# Create a table of the number of observations per subject in the longitudinal data
m = table(Long$SubjectID)
ID = unique(Long$SubjectID)
n = dim(Surv)[1]

# Initialize a matrix to store the history of preterm delivery
HH = matrix(NA, n, max(m))
for(i in 1:n){
  HH[i, 1:length(Long$History_of_preterm_delivery[Long$SubjectID == ID[i]])] = Long$History_of_preterm_delivery[Long$SubjectID == ID[i]]
}

# Add the history of preterm delivery to the survival data
Surv$History = HH[,1]

# Log-transform the library size in the longitudinal data
Long$lLibrarySize = log(Long$LibrarySize)

# Set seed for reproducibility and define training subjects
set.seed(3)
INDTRAIN = unique(Long$SubjectID)[1:30]

# Split the longitudinal and survival data into training and validation sets
dataLong_t <- subset(Long, Long$SubjectID %in% INDTRAIN)
dataSurv_t <- subset(Surv, Surv$SubjectID %in% INDTRAIN)
dataLong_v <- subset(Long, !(Long$SubjectID %in% INDTRAIN))
dataSurv_v <- subset(Surv, !(Surv$SubjectID %in% INDTRAIN))

# Fit a Zero-Inflated Negative Binomial model with shared random effects
ZNB <- ZISRE(
  FixedY = Prevotella ~ GWColl + Preeclampsia, RandomY = ~1, 
  GroupY = ~SubjectID,
  FixedZ = ~ GWColl + Preeclampsia, RandomZ = ~1, GroupZ = ~SubjectID,
  formSurv = Surv(GWDels, Deliverys) ~ Preeclampsias + non_hispanic_whites + high_incomes + History,
  IStructure = FALSE,
  obstime = "GWColl", offset = "lLibrarySize",
  dataLong = dataLong_t, dataSurv = dataSurv_t,
  n.chains = 2,
  n.iter = 100000, n.burnin = 80000, n.thin = 90, family = "NB"
)
```
The output of the package for estimation is as follows:

```
ZNB$Estimation
$Count_model
                         Est         SD         L_CI        U_CI      Rhat
(Intercept)      -4.86258618 0.62524654 -6.152726826 -3.66164495 1.0013049
GWColl            0.03177529 0.01472073  0.002814288  0.06079314 0.9999121
PreeclampsiaTRUE  1.81592976 0.92343311  0.038587236  3.63966332 1.0004539
Dispersion        0.25422084 0.03357116  0.189130248  0.32194203 1.0000550

$Zero_inflated_model
                         Est         SD       L_CI         U_CI      Rhat
(Intercept)      -1.35877745 0.69040091 -2.7533507 -0.034178881 1.0008752
GWColl           -0.02670194 0.01760999 -0.0615588  0.007203673 1.0000727
PreeclampsiaTRUE -1.06631876 1.05515926 -3.1657119  1.045817120 0.9999143

$Survival_model
                                 Est        SD         L_CI        U_CI      Rhat
(Intercept)             -71.49448770 9.8275798 -91.20222411 -51.5678938 1.2787956
PreeclampsiasTRUE         0.99378756 0.5299029  -0.06832966   2.0214544 1.0107505
non_hispanic_whitesTRUE  -0.03862752 0.5189329  -1.08372199   0.9497620 0.9999063
high_incomesTRUE          0.13843219 0.4885124  -0.80941644   1.0838178 0.9999214
History                   0.82286005 0.5126226  -0.18423381   1.8216419 1.0009330
(Intercept)               0.21072087 0.2102530  -0.16876570   0.6485199 1.0036064
(Intercept)               0.12067740 0.2008500  -0.25576124   0.5435720 1.0007091
Scale                    19.40354831 2.6676309  13.98009625  24.7394006 1.2805303

$D
$D$Est
            (Intercept) (Intercept)
(Intercept)    4.083613   -3.196732
(Intercept)   -3.196732    4.535441

$D$SD
            (Intercept) (Intercept)
(Intercept)    1.403819    1.230802
(Intercept)    1.230802    1.945547

$D$L_CI
            (Intercept) (Intercept)
(Intercept)    2.041110   -6.117572
(Intercept)   -6.117572    2.006391

$D$U_CI
            (Intercept) (Intercept)
(Intercept)    7.438937   -1.421713
(Intercept)   -1.421713    9.567356

$D$Rhat
            (Intercept) (Intercept)
(Intercept)    1.001243    1.001984
(Intercept)    1.001984    1.000450
```

Compute risk prediction, as well as AUC and BS

```
# Compute risk prediction, as well as AUC and BS
s1 = 35.5
t1 = 1

DD <- DP_SRE(
  object = ZNB, s = s1, t = t1, offset = "lLibrarySize", n.chains = 1, n.iter = 6000, n.burnin = 3000,
  n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
)

Criteria(
  s = s1, t = t1, Survt = dataSurv_v$GWDels,
  CR = dataSurv_v$Deliverys, P = DD$DP$est, cause = 1
)$Cri

s1 = 36

DD <- DP_SRE(
  object = ZNB, s = s1, t = t1, offset = "lLibrarySize", n.chains = 1, n.iter = 6000, n.burnin = 3000,
  n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
)

Criteria(
  s = s1, t = t1, Survt = dataSurv_v$GWDels,
  CR = dataSurv_v$Deliverys, P = DD$DP$est, cause = 1
)$Cri

s1 = 36.5

DD <- DP_SRE(
  object = ZNB, s = s1, t = t1, offset = "lLibrarySize", n.chains = 1, n.iter = 6000, n.burnin = 3000,
  n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
)

Criteria(
  s = s1, t = t1, Survt = dataSurv_v$GWDels,
  CR = dataSurv_v$Deliverys, P = DD$DP$est, cause = 1
)$Cri

# Compute credible intervals for the predictions
DD <- DP_SRE_CI(
  object = ZNB, s = s1, t = t1, offset = "lLibrarySize", mi = 2, n.chains = 1, n.iter = 60, n.burnin = 30,
  n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
)
```

```
# Generate dynamic predictions for patient #10055
par(mfrow=c(2,1))

DPplot2(ZNB,
        s = 30, id_new = 10055, by = 1, mi = 10, digits=0,
        Marker_lab="Prevotella", Time_lab="Time (Week)", 
        offset = "lLibrarySize",
        n.chains = 1, n.iter = 2000, n.burnin = 1000,
        dataLong = dataLong_v, dataSurv = dataSurv_v
)
title("ZINB", font.main=2, cex.main=1)

DPplot2(ZNB,
        s = 35, id_new = 10055, by = 1, mi = 10, digits=0,
        Marker_lab="Prevotella", Time_lab="Time (Week)", 
        offset = "lLibrarySize",
        n.chains = 1, n.iter = 2000, n.burnin = 1000,
        dataLong = dataLong_v, dataSurv = dataSurv_v
)
title("ZINB", font.main=2, cex.main=1)
```

Finally, the individual dynamic prediction (DP) is as follows:

<img src="/Figure/preg_plot.png" alt="Description" width="600" height="500">

The R code for ZIP and ZIGP follows the same structure as described above.











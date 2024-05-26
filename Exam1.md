R code for ZINB joint modeling
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
  n.iter = 2000, n.burnin = 1000, n.thin = 90, family = "NB"
)
```

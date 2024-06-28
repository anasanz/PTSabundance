## -------------------------------------------------
##                1. M0: Constant p in time
##                    Sex-specific p
## ------------------------------------------------- 

rm(list=ls())

library(jagsUI)
library(rjags)

## ---- Data ----

setwd("D:/MargSalas/Ganga/Data/CMR")

load("cr_sandgrouse_2022_new.RData")
capt.hist <- as.data.frame(capt.hist$ch)
colnames(capt.hist)[1] <- "ch"
n_oc <- 5

data_ganga <- matrix(data = as.numeric(do.call("rbind", strsplit(as.character(capt.hist$ch), "", fixed = TRUE))), nrow = nrow(capt.hist), ncol = n_oc)
T = ncol(data_ganga)

# Augment data set by 150 potential individuals
nz <- 150
yaug <- rbind(data_ganga, array(0, dim = c(nz, T)))

# Load sex

setwd("D:/MargSalas/Ganga/Data/CMR")
load("id_sex_sandgrouse_2022.RData") # It is in the same order than capt.hist
id_sex$sex[which(id_sex$sex == "F")] <- 0
id_sex$sex[which(id_sex$sex == "M")] <- 1
id_sex$sex[which(id_sex$sex == "X")] <- NA

sex <- as.numeric(id_sex$sex)
sexAug <- c(sex, rep(NA, nz))

# Specify model in BUGS language
setwd("D:/MargSalas/Ganga/Data/CMR")
sink("model_m0_pSex.txt")
cat("
model {

# Priors
omega ~ dunif(0, 1)
psi ~ dunif(0, 1)

for (s in 1:2){
  p[s] ~ dunif(0, 1)
}

# Likelihood
for (i in 1:M){
   sex[i] ~ dbern(psi)	
   z[i] ~ dbern(omega)			# Inclusion indicators (probability that exists)
   for (j in 1:T){
      yaug[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i] * p[sex[i]+1]		# Can only be detected if z=1
      } #j
   } #i

# Derived quantities
N <- sum(z[])
}",fill = TRUE)
sink()

# Bundle data
win.data <- list(yaug = yaug, M = nrow(yaug), T = ncol(yaug), sex = sexAug)


# Initial values
inits <- function() list(z = rep(1, nrow(yaug)), p = runif(2, 0, 1),
                         sex = c(rep(NA,nrow(data_ganga)), rbinom(nz,1,0.5)))

# Parameters monitored
params <- c("N", "p", "omega", "psi", "z", "sex")

# MCMC settings
ni <- 2500
#nt <- 2
nt <- 1
nb <- 500
nc <- 3

# Call WinBUGS from R (BRT <1 min)
out_sim_M0_psex <- jags(win.data, inits, params, "model_m0_pSex.txt", n.chains = nc, 
                        n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

## -------------------------------------------------
##                 Results
## ------------------------------------------------- 

# Estimate number of individuals of each sex

zzall <- out_sim_M0_psex$sims.list$z

zzfem <- out_sim_M0_psex$sims.list$z
zzfem[!out_sim_M0_psex$sims.list$sex %in% c(1)] <- 3 # Only females (set males as dead)

zzmal <- out_sim_M0_psex$sims.list$z
zzmal[!out_sim_M0_psex$sims.list$sex %in% c(0)] <- 3 # Only Males (set females as dead)

N <- matrix(NA, nrow = dim(out_sim_M0_psex$sims.list$z)[1], ncol = 3)

for(ite in 1:dim(out_sim_M0_psex$sims.list$z)[1]){
  
  which.alive.all <- which(zzall[ite,] == 1) # Select only the individuals alive (z=1)
  N[ite,1] <- length(which.alive.all)
  
  which.alive.fem <- which(zzfem[ite,] == 1) # Select only the individuals alive (z=1)
  N[ite,2] <- length(which.alive.fem) 
  
  which.alive.mal <- which(zzmal[ite,] == 1) # Select only the individuals alive (z=1)
  N[ite,3] <- length(which.alive.mal) 
}

(ab <- c(paste("total = ", round(mean(N[,1]),0), 
               "; female = ", round(mean(N[,2]),0), 
               "; male = ", round(mean(N[,3]),0), sep = "")))

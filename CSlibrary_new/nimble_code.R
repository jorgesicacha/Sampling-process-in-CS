library(nimble)
source("estpar.R")
load("data_for_nimble.RData")

code <-nimbleCode({
  #prior for omega
  for(i in 1:nspecies){
    alpha[i] ~ dexp(1)
    omega[i, 1:nspecies] ~ ddirch(alpha = alpha[1:nspecies])
  }
  
  r[1:nsite,1:38] <- nimble_INLA(omega[1:nspecies,1:nspecies]) #change 38 to a constant to be specified
  
  for(i in 1:nsite){
   # for(j in 1:nspecies){
    lambda_obs[i,1] <- r[i,7] + r[i,8]*true_cov[i] + r[i,13]+ 
      r[i, 9] + r[i,10]*bias_cov[i]+ r[i,14] - log(1+exp(r[i, 9] + r[i,10]*bias_cov[i]+ r[i,14])) + 
      r[i, 11] + r[i,12]*det_cov[i] - log(1+exp(r[i, 11] + r[i,12]*det_cov[i]))  
  # Second species
    lambda_obs[i,2] <- r[i,15] + r[i,16]*true_cov[i] + r[i,21]+ 
      r[i, 17] + r[i,18]*bias_cov[i]+ r[i,22] - log(1+exp(r[i, 17] + r[i,18]*bias_cov[i]+ r[i,22])) + 
      r[i, 19] + r[i,20]*det_cov[i] - log(1+exp(r[i, 19] + r[i,20]*det_cov[i]))  
    # Third species
    lambda_obs[i,3] <- r[i,23] + r[i,24]*true_cov[i] + r[i,29]+ 
      r[i, 25] + r[i,26]*bias_cov[i]+ r[i,30] - log(1+exp(r[i, 25] + r[i,26]*bias_cov[i]+ r[i,30])) + 
      r[i, 27] + r[i,28]*det_cov[i] - log(1+exp(r[i, 27] + r[i,28]*det_cov[i])) 
    # Third species
    lambda_obs[i,4] <- r[i,31] + r[i,32]*true_cov[i] + r[i,37]+ 
      r[i, 33] + r[i,34]*bias_cov[i]+ r[i,38] - log(1+exp(r[i, 33] + r[i,34]*bias_cov[i]+ r[i,38])) + 
      r[i, 35] + r[i,36]*det_cov[i] - log(1+exp(r[i,35] + r[i,36]*det_cov[i]))  
    
    # }
  }

  log(lambda[1:nsite, 1:nspecies]) <- lambda_obs[1:nsite, 1:nspecies]
  
  #Proportion of lambdas
  
  
  # Proportion for the multinomial distribution
  for(site.tag in 1:nsite){
    for(spe.tag in 1:nspecies){
      prop[site.tag,spe.tag] <- (lambda[site.tag, spe.tag])/sum(lambda[site.tag, 1:nspecies])
    }
  }
  
  
  # True data
  for(site.tag in 1:nsite){
      C[site.tag] ~ dcat(prop[site.tag,1:nspecies])
  }
  
  # Reported species
    for(site.tag in 1:nsite){
      Y[site.tag] ~ dcat(omega[C[site.tag],1:nspecies])
    }
  
})


## Parameterising the nimble model

#Data
inla_data <- list(Y=data_df$Y, 
                  C = data_df$C,
                true_cov = data_df$eco_cov,
                  bias_cov=data_df$samp_cov,
                det_cov= data_df$det_cov)

#Constants
const <- list(nspecies=length(unique(data_df$C)),
              nsite = length(data_df$C)
              )

# Initial values
idm_inits <- function(){list(beta=c(1,1)
                             
)
}

initsList <- idm_inits()

#Putting all together for the creating compilation
modelInfo <- list(
  code = code,
  constants = const,
  data = inla_data)#,
# inits = initsList
#)

#Create the model in nimble
mwtc <- nimbleModel(code, 
                    data = inla_data,
                    constants = const#, 
                    #inits = initsList
)
#library(igraph)
#plot(mwtc$modelDef$graph)

# Create the model in C
Cmwtc <- compileNimble(mwtc,showCompilerOutput = FALSE) #Have issues compiling


mcmcconf <- configureMCMC(Cmwtc, monitors = c("omega"))

Rmcmc <- buildMCMC(mcmcconf, 
                   enableWAIC =FALSE)

# Compile 
cmcmc <- compileNimble(Rmcmc, 
                       project = Cmwtc,
                       resetFunctions = TRUE)

# Run the MCMC
mcmc.out <- runMCMC(cmcmc, 
                    niter = 500,
                    nchains = 3,
                    # nburnin = 2500,
                    #inits = initsList,
                    #thin =10, 
                    setSeed = TRUE, 
                    samples=TRUE, 
                    samplesAsCodaMCMC = TRUE, 
                    summary = TRUE, 
                    WAIC = FALSE)

#Output from the MCMC
output <- mcmc.out$summary$all.chains
output

library(nimble)
library(parallel)
#source("estpar2.R")
#load("data_for_nimble.RData")
ii <- 0
listout <- list()

library(nimble)
code <-nimbleCode({
  #prior for omega
  for(i in 1:nspecies){
    #alpha[i] ~ dexp(1)
    omega[i, 1:nspecies] ~ ddirch(alpha = alpha[1:nspecies])
  }
  
  # omega[1,1] ~ dunif(0.2,1)
  # omega[2,2] ~ dunif(0.2,1)
  # 
  # omega[1,2] <- 1-omega[1,1]
  # omega[2,1] <- 1-omega[2,2]
  
  r[1:nsite,1:20] <- nimble_INLA(omega[1:nspecies,1:nspecies]) #change 38 to a constant to be specified
  
  for(i in 1:nsite){
    # for(j in 1:nspecies){
    log(lambda_obs[i,1]) <- r[i,5] + r[i,6]*true_cov[i] + r[i,11]+ 
      r[i, 7] + r[i,8]*bias_cov[i]+ r[i,12] - log(1+exp(r[i, 7] + r[i,8]*bias_cov[i]+ r[i,12])) + 
      r[i, 9] + r[i,10]*det_cov[i] - log(1+exp(r[i, 9] + r[i,10]*det_cov[i]))  
    # Second species
    log(lambda_obs[i,2]) <- r[i,13] + r[i,14]*true_cov[i] + r[i,19]+ 
      r[i, 15] + r[i,16]*bias_cov[i]+ r[i,20] - log(1+exp(r[i, 15] + r[i,16]*bias_cov[i]+ r[i,20])) + 
      r[i, 17] + r[i,18]*det_cov[i] - log(1+exp(r[i, 17] + r[i,18]*det_cov[i]))  
    
  }
  
  lambda[1:nsite, 1:nspecies] <- lambda_obs[1:nsite, 1:nspecies]
  
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
              nsite = length(data_df$C),
              alpha=c(1,1)
)

# Initial values
  idm_inits <- function(){list(omega= matrix(c(0.5, 0.5,
                                               0.5, 0.5),
                                             nrow=2, ncol=2, byrow = TRUE)
                               
  )
  }
  
  initsList <- idm_inits()
  
  #Putting all together for the creating compilation
  modelInfo <- list(
    code = code,
    constants = const,
    data = inla_data,
   inits = initsList
  )
  
  #Create the model in nimble
  mwtc <- nimbleModel(code, 
                      data = inla_data,
                      constants = const, 
                      inits = initsList
  )
  #library(igraph)
  #plot(mwtc$modelDef$graph)
  
  # Create the model in C
  Cmwtc <- compileNimble(mwtc,showCompilerOutput = FALSE) #Have issues compiling
  
  
  mcmcconf <- configureMCMC(Cmwtc, print=TRUE, useConjugacy=FALSE,monitors = c("omega","r"))
  
  mcmcconf$removeSamplers(c("omega[1,1:2]","omega[2,1:2]"))
  mcmcconf$addSampler(c("omega[1,1:2]"), "RW_dirichlet", adaptive=TRUE, scale=3)
  mcmcconf$addSampler(c("omega[2,1:2]"), "RW_dirichlet", adaptive=TRUE, scale=3)
  #mcmcconf$removeSamplers(c('omega'))
  
  Rmcmc <- buildMCMC(mcmcconf) 
                     #enableWAIC =FALSE)
  
  # Compile 
  cmcmc <- compileNimble(Rmcmc, 
                         project = Cmwtc,
                         resetFunctions = FALSE)

# Run the MCMC
library(pbapply)

  mcmc.out <- runMCMC(cmcmc, 
                            niter = 100,
                            # nburnin = 2500,
                            # inits = initsList,
                            #thin =10, 
                            #setSeed = x, 
                            samples=TRUE, 
                            samplesAsCodaMCMC = TRUE, 
                            summary = TRUE, 
                            WAIC = FALSE)

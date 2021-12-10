library(nimble)


code <-nimbleCode({
  #prior for omega
  for(i in 1:nspecies){
    alpha[i] ~ dexp(1)
    omega[i, 1:nspecies] ~ ddirch(alpha = alpha[1:nspecies])
  }
  
  r[1:24] <- nimble_INLA(omega)
  
  for(i in 1:nsite){
    for(j in 1:nspecies){
    # lambda_obs[i,j] <- r[j%%]
    }
  }

  
  inla.res[1:2] <- nimble_INLA(x_obs[1:N,1:2],y_obs[1:N],beta[1:2])
  #for(j in 1:2){
  beta[1] ~ dnorm(0,sd=1)  
  beta[2] ~ dnorm(0,sd=1)
  #}
  for(i in 1:N){
    linpred[i] <- inla.res[1]+ beta[1]*x[i,1]+ beta[2]*x[i,2]
    y[i] ~ dnorm(linpred[i],inla.res[2])
  }
  #Note that x=x_obs, y=y_obs
  #I don't know if I am doing right
  #But that is a trick I have employed
})


## Parameterising the nimble model

#Data
inla_data <- list(y_obs=df$y, 
                  x_obs = df$x,
                  y = df$y,
                  x=df$x)

#Constants
const <- list(N = 100)

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


mcmcconf <- configureMCMC(Cmwtc, monitors = c("beta", "inla.res"))

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

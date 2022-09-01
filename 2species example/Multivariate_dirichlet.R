library(nimble)
library(parallel)
source("estpar.R")
#load("data_for_nimble.RData")
ii <- 0
listout <- list()



myRW_dirichlet <- nimbleFunction(
  #name = 'sampler_RW_dirichlet',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## control list extraction
    adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
    adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
    adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
    scaleOriginal       <- extractControlElement(control, 'scale',               1)
    ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes <- model$getDependencies(target) #stochOnly = TRUE
    calcNodesNoSelf <- model$getDependencies(target,self = FALSE)
    isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)   ## should be made faster
    calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
    calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
    ## numeric value generation
    d <- length(targetAsScalar)
    thetaVec         <- rep(0, d)
    scaleVec         <- 1#scaleOriginal
    timesRan         <- 0
    timesAcceptedVec <- 0
    timesAdapted     <- 0
    optimalAR        <- 0.44
    gamma1           <- 0
    ## checks
    if(length(model$expandNodeNames(target)) > 1)    stop('RW_dirichlet sampler only applies to one target node')
    if(model$getDistribution(target) != 'ddirch')    stop('can only use RW_dirichlet sampler for dirichlet distributions')
  },
  run = function() {
    if(thetaVec[1] == 0)   thetaVec <<- values(model, target)   ## initialization
    alphaVec <- model$getParam(target, 'alpha')
    #for(i in 1:d) {
      currentValue <- thetaVec
      propLogScale <- rnorm(d, mean = 0, sd = scaleVec)
      propValue <- currentValue * exp(propLogScale)
      if(propValue[1] != 0) {
        thetaVecProp <- thetaVec
        thetaVecProp <- propValue
        values(model, target) <<- thetaVecProp / sum(thetaVecProp)
        logMHR <- alphaVec*propLogScale + currentValue - propValue + model$calculateDiff(calcNodesNoSelf)
        jump <- decide(logMHR[1])
      } else jump <- FALSE
      if(adaptive & jump)   timesAcceptedVec <<- timesAcceptedVec + 1
      if(jump) {
        thetaVec <<- thetaVecProp
        nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
        nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
      } else {
        nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
        nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
        nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
      }
      model$calculate(target)                                                             ## update target logProb
      nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProbOnly = TRUE)    ##
    #}
    if(adaptive) {
      timesRan <<- timesRan + 1
      if(timesRan %% adaptInterval == 0) {
        acceptanceRateVec <- timesAcceptedVec / timesRan
        timesAdapted <<- timesAdapted + 1
        gamma1 <<- 1/((timesAdapted + 3)^adaptFactorExponent)
        adaptFactorVec <- exp(10 * gamma1 * (acceptanceRateVec - optimalAR))
        scaleVec <<- scaleVec * adaptFactorVec
        timesRan <<- 0
        timesAcceptedVec <<- 0
      }
    }
  },
  methods = list(
    reset = function() {
      thetaVec         <<- numeric(d, 0)
      scaleVec         <<- scaleOriginal
      timesRan         <<- 0
      timesAcceptedVec <<- 0
      timesAdapted     <<- 0
      gamma1           <<- 0
    }
  )
)


assign('myRW_dirichlet', myRW_dirichlet, envir = .GlobalEnv)





dmydirch <- nimbleFunction(
  run = function(x = double(1), alpha = double(1), 
                 log = integer(0, default = 0)) {
    returnType(double(0))
    logProb <- sum(lgamma(alpha)) - lgamma(sum(alpha)) + 
      sum((alpha -1) * log(x))
    if(log) return(logProb)
    else return(exp(logProb))
  })

rmydirch <- nimbleFunction(
  run = function(n = integer(0), alpha = double(1)) {
    returnType(double(1))
    if(n != 1) print("rdirch only allows n = 1; using n = 1.")
    p <- rdirch(1, alpha)
    return(p)
  })

# omega
sum_fnx <- function(omega){
  ii  <- get("ii",envir =  parent.frame())
  ii <- assign("ii",ii+1,envir = parent.frame())
  
  print(ii)
  return(return(sum(omega)))
}

nimble_omega <- nimbleRcall(
  prototype = function(
    omega=double(1)#x is a matrix 
    # beta is a vector
  ) {},
  returnType = double(0), # outcome is a vector
  Rfun = 'sum_fnx'
)

#omega <- nimb

#my_laplace <- buildLaplace()



library(nimble)
library(dirmult)
code <-nimbleCode({
  #prior for omega
  #for(i in 1:nspecies){
  #  alpha[i] ~ dexp(1)
  # }
  
  #for(i in 1:nspecies){
  omega[1, 1:nspecies] ~ ddirch(alpha = alpha[1:nspecies])
  
  #}
  
  #p11 ~ dunif(0.2,1)
  #p22 ~ dunif(0.2,1)
  
  #omega[1:nspecies, 1:nspecies] <- nimble_omega(p11, p22)
  
  r[1,1] <- nimble_omega(omega[1,1:nspecies]) #change 38 to a constant to be specified
  
  for(site.tag in 1:nsite){
    C[site.tag] ~ dpois(r[1,1])
  }
  # Reported species
  for(site.tag in 1:nsite){
    Y[site.tag] ~ dcat(omega[1,1:nspecies])
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
              alpha=rep(1, length(unique(data_df$C)))
)

# Initial values
idm_inits <- function(){list(
omega = matrix(c(0.5, 0.5), 1,2)
)
}

initsList <- idm_inits()

#Putting all together for the creating compilation
modelInfo <- list(
  code = code,
  constants = const,
  data = inla_data#,
  #inits = initsList
)

#Create the model in nimble
mwtc <- nimbleModel(code, 
                    data = inla_data,
                    constants = const, 
                    inits = initsList)
#)
#library(igraph)
#plot(mwtc$modelDef$graph)

# Create the model in C
Cmwtc <- compileNimble(mwtc,
                       showCompilerOutput = FALSE) #Have issues compiling


mcmcconf <- configureMCMC(Cmwtc, 
                          print=TRUE, 
                          useConjugacy=FALSE,
                          monitors = c("omega","r"))

mcmcconf$removeSamplers(c("omega[1,1:2]"))
mcmcconf$addSampler(c("omega[1,1:2]"), "myRW_dirichlet")
#mcmcconf$addSampler(c("omega[2,1:2]"), "myRW_dirichlet")




Rmcmc <- buildMCMC(mcmcconf) 
#enableWAIC =FALSE)

# Compile 
cmcmc <- compileNimble(Rmcmc, 
                       project = Cmwtc,
                       resetFunctions = FALSE)

# Run the MCMC
#library(pbapply)

mcmc.out <- runMCMC(cmcmc, 
                    niter = 20,
                    #nburnin = 5000,
                    inits = initsList,
                    # thin =100, 
                    #setSeed = x, 
                    samples=TRUE, 
                    samplesAsCodaMCMC = TRUE, 
                    summary = TRUE, 
                    WAIC = FALSE)
#save(mcmc.out, file="estimates.RData")

library(nimble)
est_par <- function(omega){
  
  if(any(omega) < 0 | any(omega) > 1) stop("Components of omega should be positive")
  
  
  #inla.setOption(pardiso.license= "pardiso.lic")
  #control.compute=list(openmp.strategy="pardiso.parallel")
  #inla.pardiso.check()
  # x<- get("x")
  # csdata <- get("csdata",envir = eval(get(paste0("a",x))))
  # cssampdata <- get("cssampdata",envir = eval(get(paste0("a",x))))
  # detdata <- get("detdata",envir = eval(get(paste0("a",x))))
  # covs <- get("covs",envir = eval(get(paste0("a",x))))
  # region <- get("region",envir = eval(get(paste0("a",x))))
  # mesh <- get("mesh", envir = eval(get(paste0("a",x))))
  # listout <- get("listout",envir = eval(get(paste0("a",x))))
  # ii  <- get("ii",envir =  eval(get(paste0("a",x))))
  # ii <- assign("ii",ii+1,envir = eval(get(paste0("a",x))))
  # 
  # print(ii)
  # print(x)
  # 
  
  csdata <- get("csdata",envir = parent.frame())
  cssampdata <- get("cssampdata",envir = parent.frame())
  detdata <- get("detdata",envir = parent.frame())
  covs <- get("covs",envir = parent.frame())
  region <- get("region",envir = parent.frame())
  mesh <- get("mesh", envir = parent.frame())
  listout <- get("listout",envir = parent.frame())
  ii  <- get("ii",envir =  parent.frame())
  ii <- assign("ii",ii+1,envir = parent.frame())
  
   print(ii)
  # print(x)
  # Organizing the inputs
  nspecies = nrow(omega)
  tmp <- csdata$classification
  Eco_PPFinal_detect <- list()
  for(i in 1:nspecies){Eco_PPFinal_detect[[i]] <- tmp[which(tmp$error==i),]}
  Samp_PPFinal <- cssampdata
  data_det_spframe <- detdata
  # cov1.spix <- covs[[1]]
  # cov2.spix <- covs[[2]]
  # cov3.spix <- covs[[3]]
  aa <- region
  # mesh <- mesh
  # spdes <- spdes$spdes
  # spde2 <- spdes$spde2
  
  
  #Defining components of the model
  
  
  cmp1 <- list()
  for(i in 1:nspecies){
    cmp1[[i]] <- (paste0("+ beta0",i,"(1)", 
                         "+ beta0thin(1)",
                         "+ beta0det" ,i, "(1)",
                         "+w1",i,"(main = coordinates, model =","spdes[[",i, "]])", 
                         "+ w2(main = coordinates, model = spde2)+",
                         "cov1",i, "(main=cov1.spix,model='linear') +",
                         "cov2(main=cov2.spix,model='linear')+",
                         "cov3",i, "(main=cov3.spix,model='linear')"))
  }
  cmp <- as.formula(paste0("~ -1",do.call("paste0",cmp1)))
  
  fun <- function(x,y,z){
    -log(1+exp((x+y+z)))
  }
  
  fun1 <- function(x,y){
    -log(1+exp((x+y)))
  }
  
  
  
  fun21 <-  function(a,b,c,d,e,f){
    ret <- log(omega[1,1]*plogis(a+b+c) + omega[2,1]*plogis(d+e+f))
    #ret <- log((exp(a+b+c)/(1+exp(a+b+c)))*omega1)
    return(ret)
  }
  
  fun22 <-  function(a,b,c,d,e,f){
    ret <- log(omega[1,2]*plogis(a+b+c) + omega[2,2]*plogis(d+e+f))
    return(ret)
  }
  

  
  
  lik1 <- lik2 <- lik3 <- list()
  
  for(i in 1:nspecies){
    lik1[[i]] <- inlabru::like("cp",
                               formula = as.formula(paste0("coordinates ~ beta0",i,"  + cov1",i," + w1",i, "+ beta0thin + cov2 + w2 + fun(beta0thin,cov2,w2)+ beta0det",i,"+cov3",i,"+fun1(beta0det",i,", cov3",i,")",
                                                           "+fun2",i,"(beta01, cov11, w11, beta02, cov12, w12)")),
                               data = Eco_PPFinal_detect[[i]],
                               #components = cmp,
                               domain = list(coordinates = mesh),
                               samplers = aa)
    lik2[[i]] <- inlabru::like("cp",
                               formula = coordinates ~ beta0thin + cov2 + w2,
                               data = Samp_PPFinal,
                               #components = cmp,
                               domain = list(coordinates = mesh),
                               samplers = aa)
    lik3[[i]] <- inlabru::like("binomial",
                               formula = as.formula(paste0("detdata",i," ~ beta0det",i," + cov3",i)),
                               data = data_det_spframe[[i]],
                               #components = cmp,
                               domain = list(coordinates = mesh),
                               samplers = aa)
  }
  
  
  #str(lik1[[1]])
  
  
  # 
  # 
  
  predpoints <- expand.grid(x=seq(0,3,length.out = 128),y=seq(0,3,length.out = 128))
  cov1.pred <- cos(predpoints$x) - sin(predpoints$y - 2)
  cov2.pred <- cos(2*predpoints$x) - sin(2*predpoints$y-4)
  cov3.pred <- (predpoints$x/2)^2+(predpoints$y/2)^2
  
  
  # par(mfrow=c(1,3))
  # plot(pred.rast.median.bru)#,zlim=c(0.9,4.1))
  # plot(pred.rast.sd.bru)
  # plot(w1.rast)
  
  inlabru:::iinla.setOption("iinla.verbose", TRUE)
  #fit2 <- list()
  pred.median.eco <- list()
  pred.sd.eco <- list()
  pred.median.samp <- list()
  pred.sd.samp <- list()
  pred.median.det <- list()
  pred.sd.det <- list()
  
  cov1.spix$cov11 <- cov1.spix$cov12 <- cov1.spix$layer
  cov3.spix$cov31 <- cov3.spix$cov32  <- cov3.spix$layer
  #for(i in 1:nspecies){
  # names(cov1.spix)<- paste0("cov1",i)
  #nemes(cov1.spix) <- "cov11"
  names(cov2.spix)<- "cov2"
  #names(cov3.spix) <- paste0("cov3",i)
  fit2 <- inlabru::bru(cmp, lik1[[1]], lik1[[2]],
                       lik2[[1]],
                       lik3[[1]],lik3[[2]],
                       options = list(control.inla = list(strategy = "gaussian",
                                                          int.strategy = "eb"),
                                      bru_method = list(
                                        taylor = "pandemic",
                                        search = "all",
                                        factor = (1 + sqrt(5)) / 2,
                                        rel_tol = 0.99,
                                        max_step = 2,
                                        lin_opt_method = "onestep"
                                      ), #change to 0.01
                                      bru_max_iter=5)) #change to 50
  
  #postsamples <- INLA::inla.posterior.sample.eval(fun=function(x){x},samples=INLA::inla.posterior.sample(n = 1000,result = fit2,intern = FALSE))
  
  if(length(listout)==0){assign("listout",fit2)}
  else{assign("listout",rbind(listout,fit2))}
  
  #assign("listout",rbind(listout,fit2),envir=parent.frame())
  
  plots = FALSE
  if(plots == TRUE){
    require(ggpubr)
true_values <- c(0.8, 1.3, 2, 1.5, -1.5, -2, 2.5, -0.3, -0.12, -0.5)
title_names <- noquote(c("beta[01]", "alpha[01]", "gamma[01]", 
                 "beta[11]", "alpha[11]", "gamma[11]",
                 "beta[02]", "gamma[02]",
                 "beta[12]",  "gamma[12]"
                 ))
    pp <- lapply(as.list(c(seq(1:length(fit2$names.fixed)))), function(x){
      plot(fit2, fit2$names.fixed[[x]])+
        geom_vline(xintercept = true_values[x], linetype = 2, col = "red")+
        theme_bw()+
        labs(title = parse(text = title_names[x]))+
        xlab("")+
        ylab("")
      
      #p + ggtitle(label =expression(title_names[x]) )
      #p + annotate("text",x = 2, y = 0.3, 
        #           parse = TRUE, label = title_names[[x]])
      })
    
    ggarrange(pp[[1]], pp[[2]], pp[[3]], pp[[4]], pp[[5]],
              pp[[6]], pp[[7]], pp[[8]], pp[[9]], pp[[10]])
    
    names_random <- c("Range for w11", "Stdev for w11",
                              "Range for w2", "Stdev for w2",
                              "Range for w12", "Stdev for w12")
    true_values_r <- c(1.2, 0.2, 2.5, 0.2, 2.5, 1.2 )
    ppr <- lapply(as.list(c(seq(1:length(names_random)))), function(x){
      plot(fit2, names_random[x])+
        geom_vline(xintercept = true_values_r[x], linetype = 2, col = "red")+
        theme_bw()+
        labs(title = names_random[x])+
        xlab("")+
        ylab("")
    })
    
    ggarrange(pp[[1]], pp[[2]], pp[[3]], pp[[4]], pp[[5]],
              pp[[6]], pp[[7]], pp[[8]], pp[[9]], pp[[10]],
              ppr[[1]],ppr[[2]], ppr[[3]], ppr[[4]], ppr[[5]],
              ppr[[6]] )
    
  }
  

  
  
  
  alpha0 <- alpha1 <- beta0 <- beta1 <- gamma0 <- gamma1<- c()
  tmp <- fit2$marginals.fixed
  alpha0 <- c(alpha0, INLA::inla.emarginal(function(x) x,tmp[paste0("beta0thin")][[1]]))
  alpha1 <- c(alpha1, INLA::inla.emarginal(function(x) x,tmp[paste0("cov2")][[1]]))  
  for(i in 1:nspecies){
    beta0 <- c(beta0,INLA::inla.emarginal(function(x) x,tmp[paste0("beta0",i)][[1]]))
    beta1 <- c(beta1,INLA::inla.emarginal(function(x) x,tmp[paste0("cov1",i)][[1]]))
    gamma0 <- c(gamma0,INLA::inla.emarginal(function(x) x,tmp[paste0("beta0det",i)][[1]]))
    gamma1 <- c(gamma1,INLA::inla.emarginal(function(x) x,tmp[paste0("cov3",i)][[1]]))
    #beta0 <- c(beta0,INLA::inla.emarginal(function(x) x,tmp[paste0("beta0",i)][[1]]))
    #beta0 <- c(beta0,INLA::inla.emarginal(function(x) x,tmp[paste0("beta0",i)][[1]]))
  }
  #ret <- c(beta0, beta1, alpha0, alpha1, gamma0, gamma1)
  #return(ret)
  
  ## The locations where information is needed 
  
  w11 <- sapply(fit2$marginals.random$w11,function(x){INLA::inla.emarginal(function(x) x,x)})
  w12 <- sapply(fit2$marginals.random$w12,function(x){INLA::inla.emarginal(function(x) x,x)})
  w2 <- sapply(fit2$marginals.random$w2,function(x){INLA::inla.emarginal(function(x) x,x)})
  
  
  
  obs.df <- list()
  for(i in 1:nspecies){
    obs.df[[i]] <- data.frame(coordx = lik1[[i]]$data@coords[,1],coordy = lik1[[i]]$data@coords[,2],
                              y=lik1[[i]]$data$BRU_aggregate
    )
    
    obs.df[[i]] <- obs.df[[i]][which(obs.df[[i]]$y==1),]
    names(obs.df[[i]])[3:ncol(obs.df[[i]])] <- paste0(names(obs.df[[i]])[3:ncol(obs.df[[i]])],i)
  }
  
  full_reports <- Reduce(
    function(x, y, ...) merge(x, y, all = TRUE, ...),
    obs.df
  )
  
  full_reports[is.na(full_reports)] <- 0
  
  
  A.mat <- INLA::inla.spde.make.A(mesh = mesh, loc = as.matrix(full_reports[,1:2]))
  
  w11 <- (A.mat %*% w11)[,1]
  w12 <- (A.mat %*% w12)[,1]
  w2 <- (A.mat %*% w2)[,1]
  
  pars.df <- list()
  for(i in 1:nspecies){
    pars.df[[i]] <- data.frame(coordx = full_reports$coordx,coordy = full_reports$coordy,
                               beta0 = rep(beta0[[i]],nrow(full_reports)),beta1 = rep(beta1[[i]],nrow(full_reports)),
                               alpha0 = rep(alpha0,nrow(full_reports)),alpha1 = rep(alpha1,nrow(full_reports)),
                               gamma0 = rep(gamma0[[i]],nrow(full_reports)),gamma1 = rep(gamma1[[i]],nrow(full_reports)),
                               w1=get(paste0("w1",i)),w2=w2)
    
    names(pars.df[[i]])[3:ncol(pars.df[[i]])] <- paste0(names(pars.df[[i]])[3:ncol(pars.df[[i]])],i)
  }
  
  full_pars <- Reduce(
    function(x, y, ...) merge(x, y, all = TRUE, ...),
    pars.df
  )
  
  final_df <- merge(full_reports,full_pars,by=c("coordx","coordy"),all=TRUE)
  
  return(as.matrix(final_df))                                                  
}

#x <- 0
#a0 <- as.environment(as.list(parent.frame(), all.names=TRUE))
nimble_INLA <- nimbleRcall(
  prototype = function(
    omega=double(2) #x is a matrix 
    # beta is a vector
  ) {},
  returnType = double(2), # outcome is a vector
  Rfun = 'est_par'
)

  #CnimbleINLA <- compileNimble(nimble_INLA)
  
  #Testing the compiled function. 
  #Should give the same results as fit.inla
  
  #ii <- 0
  listout <- list()
  
   class_prob <- matrix(c(0.17, 0.83,
                          0.96, 0.04),
                      nrow=2, ncol=2, byrow = TRUE)
   CnimbleINLA(class_prob)

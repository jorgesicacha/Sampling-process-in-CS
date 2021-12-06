est_par <- function(omega){
  
  csdata <- get("csdata",envir = parent.frame())
  cssampdata <- get("cssampdata",envir = parent.frame())
  detdata <- get("detdata",envir = parent.frame())
  covs <- get("covs",envir = parent.frame())
  region <- get("region",envir = parent.frame())
  mesh <- get("mesh",envir = parent.frame())
   
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
                         "+w1",i,"(map = coordinates, model =","spdes[[",i, "]])", 
                         "+ w2(map = coordinates, model = spde2)+",
                         "cov1",i, "(map=cov1.spix,model='linear') +",
                         "cov2(map=cov2.spix,model='linear')+",
                         "cov3",i, "(map=cov3.spix,model='linear')"))
  }
  cmp <- as.formula(paste0("~ -1",do.call("paste0",cmp1)))
  
  fun <- function(x,y,z){
    -log(1+exp((x+y+z)))
  }
  
  fun1 <- function(x,y){
    -log(1+exp((x+y)))
  }
  
  
  
  fun21 <-  function(a,b,c,d,e,f,g,h,i,j,k,l){
    ret <- log(omega[1,1]*plogis(a+b+c) + omega[2,1]*plogis(d+e+f) + 
                 omega[3,1]*plogis(g+h+i) + omega[4,1]*plogis(j+k+l))
    #ret <- log((exp(a+b+c)/(1+exp(a+b+c)))*omega1)
    return(ret)
  }
  
  fun22 <-  function(a,b,c,d,e,f,g,h,i,j,k,l){
    ret <- log(omega[1,2]*plogis(a+b+c) + omega[2,2]*plogis(d+e+f) + 
                 omega[3,2]*plogis(g+h+i) + omega[4,2]*plogis(j+k+l))
    return(ret)
  }
  
  fun23 <-  function(a,b,c,d,e,f,g,h,i,j,k,l){
    ret <- log(omega[1,3]*plogis(a+b+c) + omega[2,3]*plogis(d+e+f) + 
                 omega[3,3]*plogis(g+h+i) + omega[4,3]*plogis(j+k+l))
    #ret <- log((exp(a+b+c)/(1+exp(a+b+c)))*omega1)
    return(ret)
  }
  
  fun24 <-  function(a,b,c,d,e,f,g,h,i,j,k,l){
    ret <- log(omega[1,4]*plogis(a+b+c) + omega[2,4]*plogis(d+e+f) + 
                 omega[3,4]*plogis(g+h+i) + omega[4,4]*plogis(j+k+l))
    #ret <- log((exp(a+b+c)/(1+exp(a+b+c)))*omega1)
    return(ret)
  }
  
  
  lik1 <- lik2 <- lik3 <- list()
  
  for(i in 1:nspecies){
    lik1[[i]] <- inlabru::like("cp",
                      formula = as.formula(paste0("coordinates ~ beta0",i,"  + cov1",i," + w1",i, "+ beta0thin + cov2 + w2 + fun(beta0thin,cov2,w2)+ beta0det",i,"+cov3",i,"+fun1(beta0det",i,", cov3",i,")",
                                                  "+fun2",i,"(beta01, cov11, w11, beta02, cov12, w12,beta03, cov13, w13,beta04, cov14, w14)")),
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
  
  cov1.spix$cov11 <- cov1.spix$cov12 <- cov1.spix$cov13 <- cov1.spix$cov14 <- cov1.spix$layer
  cov3.spix$cov31 <- cov3.spix$cov32 <- cov3.spix$cov33 <- cov3.spix$cov34 <- cov3.spix$layer
  #for(i in 1:nspecies){
  # names(cov1.spix)<- paste0("cov1",i)
  #nemes(cov1.spix) <- "cov11"
  names(cov2.spix)<- "cov2"
  #names(cov3.spix) <- paste0("cov3",i)
  fit2 <- inlabru::bru(cmp, lik1[[1]], lik1[[2]],lik1[[3]],lik1[[4]],
              lik2[[1]],
              lik3[[1]],lik3[[2]],lik3[[3]],lik3[[4]],
              options = list(control.inla = list(strategy = "gaussian",
                                                 int.strategy = "eb"),
                             max.iter=50))
  return(fit2)
}

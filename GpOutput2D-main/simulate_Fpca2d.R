simulate.km_Fpca2d <-function(object,nsim,newdata, nugget = 0,...){
  
  # number of principal component
  nPC <- length(object)
  
  
  # number of prediction
  n<-nrow(newdata)
  
  # fpca object
  fpca <-attr(object,"fpca")
  
  # return
  res <-list()
  
  # prediction on each principal component
  p_PC <-lapply(1:nPC, function(i){
    pi <-simulate(object=object[[i]], nsim = nsim,newdata=newdata,type=type,
                  compute = compute, cond = TRUE, nugget.sim = nugget,...)
    return(pi)
  })# end p.PC
  
  
  #%%%%%%%%%%%%%
  # prediction
  #%%%%%%%%%%%%%# scores
  res = list()
  for(i in 1:nsim){
    p_pc_i = p_PC[[1]][i,]
    if(nPC > 1){for(k in 2:nPC){p_pc_i = cbind(p_pc_i, p_PC[[k]][i,])}}
    res[[i]] <- inverse_Fpca2d(cbind(p_pc_i),fpca)
  }
  
  # maps
  
  return(res)
  
}

simulate_cheap = function(model, nsim, newdata, pilot_points, covinverse = NULL, covmat1mat2 = NULL, nugget_inversion = 10^-10, nugget.sim, Z = NULL, checkNames = FALSE){
  nb_pilot = nrow(pilot_points)
  set.seed(1234)
  if(is.null(Z)){
    Z = matrix(t(simulate(model, nsim = nsim, newdata = pilot_points, nugget.sim = nugget.sim, cond = TRUE, checkNames = checkNames)), ncol = nsim)
  }
  if(is.null(covinverse)){covinverse = solve(covMatrix(model@covariance, pilot_points)$C + nugget_inversion*diag(nb_pilot))}
  if(is.null(covmat1mat2)){covmat1mat2 = covMat1Mat2(model@covariance, newdata, pilot_points)}
  
  a = matrix(rep(model@trend.coef, nsim*nrow(newdata)), ncol = nsim)
  b = matrix(covinverse%*%t(covmat1mat2), nrow = nb_pilot)
  #b = matrix(solve(K + 10^(-10)*diag(nrow(K))) %*% t(k), nrow = nb_pilot)
  #print(dim(Z))
  Z_tilde = a + t(b)%*%(Z-matrix(rep(model@trend.coef, nsim*nb_pilot), ncol = nsim))
  return(t(Z_tilde))
}

simulate_cheap.km_Fpca2d <-function(object,pilot_points, nsim,newdata, covinverse = NULL, covmat1mat2 = NULL, nugget.sim = 10^-10,nugget_inversion = 10^-10, Z = NULL, checkNames = FALSE){
  
  # number of principal component
  nPC <- length(object)
  
  
  # number of prediction
  n<-nrow(newdata)
  
  # fpca object
  fpca <-attr(object,"fpca")
  
  # return
  res <-list()
  
  # prediction on each principal component
  p_PC <-lapply(1:nPC, function(i){
    pi <-simulate_cheap(model=object[[i]], pilot_points = pilot_points, nsim = nsim,newdata=newdata, covinverse = covinverse[[i]], covmat1mat2 = covmat1mat2[[i]], nugget.sim = nugget.sim, nugget_inversion = nugget_inversion, Z = Z[[i]], checkNames = checkNames)
    return(pi)
  })# end p.PC
  
  
  res = list()
  for(i in 1:nsim){
    p_pc_i = p_PC[[1]][i,]
    if(nPC > 1){for(k in 2:nPC){p_pc_i = cbind(p_pc_i, p_PC[[k]][i,])}}
    res[[i]] <- inverse_Fpca2d(cbind(p_pc_i),fpca)
  }
  
  # maps
  
  return(res)
  
}

simulate_cheap_partial.km_Fpca2d <-function(object,pilot_points, nsim,newdata, covinverse = NULL, covmat1mat2 = NULL, nugget.sim = 10^-10,nugget_inversion = 10^-10, Z = NULL, checkNames = FALSE){
  
  # number of principal component
  nPC <- length(object)
  
  
  # number of prediction
  n<-nrow(newdata)
  
  # fpca object
  fpca <-attr(object,"fpca")
  
  # return
  res <-list()
  
  # prediction on each principal component
  p_PC <-lapply(1:nPC, function(i){
    pi <-simulate_cheap(model=object[[i]], pilot_points = pilot_points, nsim = nsim,newdata=newdata, covinverse = covinverse[[i]], covmat1mat2 = covmat1mat2[[i]], nugget.sim = nugget.sim, nugget_inversion = nugget_inversion, Z = Z[[i]], checkNames = checkNames)
    return(pi)
  })# end p.PC
  
  
  # maps
  
  return(p_PC)
  
}


simulate_partial.km_Fpca2d <-function(object,nsim,newdata, nugget = 0,...){
  
  # number of principal component
  nPC <- length(object)
  
  
  # number of prediction
  n<-nrow(newdata)
  
  # fpca object
  fpca <-attr(object,"fpca")
  
  # return
  res <-list()
  
  # prediction on each principal component
  p_PC <-lapply(1:nPC, function(i){
    pi <-simulate(object=object[[i]], nsim = nsim,newdata=newdata,type=type,
                  compute = compute, cond = TRUE, nugget.sim = nugget,...)
    return(pi)
  })# end p.PC
  
  
  #%%%%%%%%%%%%%
  # prediction

  
  return(p_PC)
  
}

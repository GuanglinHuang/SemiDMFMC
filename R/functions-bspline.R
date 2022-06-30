bs <- splines::bs

g.inverse = function(X,tol = 10^-5){
  evd = eigen(X)
  
  ev = evd$values
  evt = evd$vectors
  
  id = which(ev > tol)
  
  Mat.inv = evt[,id]%*%diag(ev[id]^-1)%*%t(evt[,id])
  
  return(Mat.inv)
}


GLR.test = function(yy,xx,mm = 1,KK = 10,rep = 100,trim=0.025,...){
  
  yy = as.matrix(yy)
  xx = as.matrix(xx)
  
  n = NCOL(yy)
  t = NROW(yy)
  
  lmInitial = lm(yy ~ xx)
  
  #(1)BSPLINE-ECSNP estimation
  
  result_bsp = coef.bspline(yy,xx,mm = mm,alpha_inl = NULL,K = KK,trim = trim,itermax = 1000)
  
  alpha_est = result_bsp$alpha_est
  e_np = result_bsp$eps_est
  e_lm = lmInitial$residuals
  
  beta_lm = t(lmInitial$coefficients[-1,])
  u_lm = lmInitial$coefficients[1,]
  
  lnL_np = 0
  lnL_lm = 0
  
  lnL_np_i = vector()
  lnL_lm_i = vector()
  
  for (nn in 1:NCOL(e_np)) {
    
    z_np = e_np[-c(1:mm),nn]/sd(e_np[-c(1:mm),nn])
    z_lm = e_lm[-c(1:mm),nn]/sd(e_lm[-c(1:mm),nn])
    
    theta_con_inl = runif(2,-1,1)
    
    constant_np = try(optim(par = theta_con_inl, function(theta){
      con = likelihood_tgc_con(z_np,theta)
      return(-con[[1]])
    },hessian = F),silent = T)
    
    constant_lm = try(optim(par = theta_con_inl, function(theta){
      con = likelihood_tgc_con(z_lm,theta)
      return(-con[[1]])
    },hessian = F),silent = T)
    
    lnL_np = lnL_np + constant_np$value - length(z_np)*log(sd(e_np[-c(1:mm),nn]))
    
    lnL_lm = lnL_lm + constant_lm$value - length(z_np)*log(sd(e_lm[-c(1:mm),nn]))
    
    lnL_np_i[nn] = constant_np$value - length(z_np)*log(sd(e_np[-c(1:mm),nn]))
    
    lnL_lm_i[nn] = constant_lm$value - length(z_np)*log(sd(e_lm[-c(1:mm),nn]))
  }
  
  GLR = (lnL_np - lnL_lm)
  
  GLR_i = (lnL_np_i - lnL_lm_i)
  
  ####bootstrap
  GLR_bs_tol = vector()
  GLR_bs_tol_i = matrix(0,rep,n)
  
  for (mmm in 1:rep) {
    
    yy_bs = matrix(NA,t,n)
    
    for (nn in 1:NCOL(e_np)) {
      
      e_star_i = as.vector(quantile(e_np[-c(1:mm),nn],runif(t)))  
      
      yy_bs[,nn] <- u_lm[nn]*rep(1,t) + xx%*%beta_lm[nn,] + e_star_i
    }
    
    
    lm_bs = lm(yy_bs ~ xx)
    
    #(1)FACE-ECSNP estimation
    
    result_bsp_bs = beta.bspline(yy_bs,xx,mm = mm,alpha = alpha_est,K = KK,trim = trim)
    
    e_np_bs = result_bsp_bs$eps_est
    e_lm_bs = lm_bs$residuals
    
    lnL_np_bs = 0
    lnL_lm_bs = 0
    lnL_np_bs_i = vector()
    lnL_lm_bs_i = vector()
    for (nn in 1:n) {
      
      z_np_bs = e_np_bs[-c(1:mm),nn]/sd(e_np_bs[-c(1:mm),nn])
      z_lm_bs = e_lm_bs[-c(1:mm),nn]/sd(e_lm_bs[-c(1:mm),nn])
      
      theta_con_inl = runif(2,-1,1)
      
      constant_np = try(optim(par = theta_con_inl, function(theta){
        con = likelihood_tgc_con(z_np_bs,theta)
        return(-con[[1]])
      },hessian = F),silent = T)
      
      constant_lm = try(optim(par = theta_con_inl, function(theta){
        con = likelihood_tgc_con(z_lm_bs,theta)
        return(-con[[1]])
      },hessian = F),silent = T)
      
      lnL_np_bs = lnL_np_bs + constant_np$value - length(z_np_bs)*log(sd(e_np_bs[-c(1:mm),nn]))
      lnL_lm_bs = lnL_lm_bs + constant_lm$value - length(z_np_bs)*log(sd(e_lm_bs[-c(1:mm),nn]))
      
      lnL_np_bs_i[nn] = constant_np$value - length(z_np_bs)*log(sd(e_np_bs[-c(1:mm),nn]))
      lnL_lm_bs_i[nn] = constant_lm$value - length(z_np_bs)*log(sd(e_lm_bs[-c(1:mm),nn]))
      
    }
    
    GLR_bs = (lnL_np_bs - lnL_lm_bs)
    
    GLR_bs_i = (lnL_np_bs_i - lnL_lm_bs_i)
    
    GLR_bs_tol = c(GLR_bs_tol,GLR_bs)
    
    GLR_bs_tol_i[mmm,] <-  GLR_bs_i
    
    print(paste0("GLR test-","Bootstrap: ",mmm))
    
  }
  
  pvalue   = 1 - length(which(GLR > GLR_bs_tol))/rep
  pvalue_i = vector()
  for (jj in 1:n) {
    pvalue_i[jj] <- 1 - length(which(GLR_i[jj] > GLR_bs_tol_i[,jj]))/rep
  }
  
  
  return(list(pvalue        = pvalue,
              pvalue_i      = pvalue_i,
              GLR           = GLR,
              GLR_i         = GLR_i,
              Boots.PDF     = GLR_bs_tol,
              Boots.PDF_i   = GLR_bs_tol_i,
              result_bsp    = result_bsp,
              result_lm     = lmInitial
  ))
  
}

calc_grid_points = function(vv,beta,num_gp)
  ## Make sure that vv only has elements 1 to n!
  ## We don't want any future values messing things up!
{
  beta = as.matrix(beta)
  zzz_temp          = vv%*%beta  
  min_zzz_temp      = min(zzz_temp) 
  max_zzz_temp      = max(zzz_temp)
  length_zzz_temp   = max_zzz_temp-min_zzz_temp;
  interval_zzz_temp = length_zzz_temp/ (num_gp-1.0)
  
  zzz_temp = matrix(NA,nrow=num_gp,ncol=1)
  for(ii in 0:(num_gp-1))
  {
    zzz_temp[ii+1] =  min_zzz_temp + ii*interval_zzz_temp;
  }
  return(as.numeric(zzz_temp))
}



interpolate_g = function(uu, theta,zz)
  ## This function estimates the intercept vector g via interpolation
  ##
  ## theta    is a (2d+2)x(num_gp)xp array of coefficient function estimates
  ## zz       are num_gp grid points that theta is estimated on
  ## uu       is a point (or vector of points) for which we would like to
  ##          interpolate over. (Formerly called z_est)
{
  d = (dim(theta)[1]-2)/2;  p = dim(theta)[3]
  g_ret =array(0, dim=c(length(uu), p, 1))
  for(k in 1:p){g_ret[,k,1] = as.matrix(approx(zz,theta[1,,k], uu, rule=2)$y)}
  
  return(g_ret)
}


interpolate_A = function(uu, theta,zz)
{
  ## This function estimates the factor loading matrix A via interpolation
  ##
  ## theta    is a (2d+2)x(num_gp)xp array of coefficient function estimates
  ## zz       are num_gp grid points that theta is estimated on
  ## uu       is the index (or vector of indicies) for which we would like to
  ##          interpolate over. (Formerly called z_est)
  
  d = (dim(theta)[1]-2)/2
  p = dim(theta)[3]
  A_ret = array(0, dim=c(length(uu), p, d))
  for(k in 1:p){for(j in 1:d){A_ret[,k,j] = approx(zz,theta[j+1,,k], uu, rule=2)$y}}     
  
  return(A_ret)
}








bs_new = function(ss,K,mm = 1){
  
  knots = seq(from = 1.00*min(ss[-c(1:mm)]),to = 1.00*max(ss[-c(1:mm)]),length.out = K-3)
  
  return(bs(ss,degree = K,intercept = F))
  
}

alpha_est_inl = function(yy,xx,alpha_inl,mm = 1,kk = 2,eps = 10^-5,itermax = 500){
  
  nn = NROW(yy)
  dd = NCOL(xx)
  pp = NCOL(yy)
  
  alpha = alpha_inl
  
  vv = matrix(NA,nrow=nn,ncol= dd)
  
  if(dd == 1){
    for(i in (mm+1):nn){if(mm==1){vv[i,] = xx[(i-1),]}else{vv[i,] = sum(xx[(i-1):(i-mm),])/mm}}
  }else{
    for(i in (mm+1):nn){if(mm==1){vv[i,] = xx[(i-1),]}else{vv[i,] = colSums(xx[(i-1):(i-mm),])/mm}}}
  
  alpha_new = rep(0,dd)  
  iter = 0
  eps = eps
  obj = 0
  while ( sum(abs(alpha_new - alpha)) > eps  ) {
    
    alpha_new = alpha
    
    ss = vv%*%alpha_new
    
    ss_tol = matrix(NA,nn,kk+1)
    
    for (i in 1:(kk+1)) {
      ss_tol[,i] = ss^(i-1)
    }
    
    xx_tol = vector()
    for (i in 1:dd) {
      xx_tol = cbind(xx_tol, ss_tol*(xx[,i]%*%t(rep(1,kk+1))))
    }
    
    reg_con = lm(yy ~ xx_tol)
    
    gamma = t(reg_con$coefficients[-1,])
    
    xx_list = list()
    
    for (i in 1:(kk+1)) {
      xx_list[[i]] = xx%*%t(gamma[,seq(1,NCOL(xx_tol),by = kk+1) + i - 1])
    }
    
    Loss_alpha = function(alpha,xx_list,yy,vv){
      xx_mat = matrix(0,nn,pp)
      for (i in 1:(kk+1)){
        xx_mat = xx_mat + matrix((vv%*%alpha)^(i-1),nn,pp,byrow = F)*xx_list[[i]]
        
      }
      
      obj = sqrt(sum(na.omit((yy - xx_mat)^2)))
      
    }
    
    alpha = rescale(optim(par = alpha,fn = function(z){Loss_alpha(z,xx_list,yy,vv)},method = "L-BFGS-B", lower = c(0,rep(-Inf,dd-1)))$par)
    
    if(alpha[1] < 0){
      alpha = - alpha
    }
    
    iter = iter + 1
    
    obj = c(obj,Loss_alpha(alpha_new,xx_list,yy,vv))
    
    if(iter > itermax){
      break
    }
    
  }
  
  
  
  return(list(alpha = alpha, obj = obj, iter = iter))
}


coef.bspline = function(yy,xx,zz = NULL,alpha_inl = NULL,mm = 1,K = 10,Kmin = 5, Kmax = 15, eps = 10^-4, itermax = 500, trim = 0.05,
                        Maxisbest = F,n.sim = 100,n.restarts = 10,LB = 0,UB = 1,trace = 1,...){
  
  LnAlpha_inl = function(alpha,vv,yy,xx,K,mm,...){
    
    alpha = rescale(alpha)
    
    pp = NCOL(yy)
    dd = NCOL(vv)
    qq = NCOL(xx)
    nn = NROW(yy)
    
    ss = vv%*%alpha
    
    bsp = bs_new(ss,K,mm)
    
    ZZ = bsp
    
    for (ii in 1:qq) {
      ZZ = cbind(ZZ,matrix(xx[,ii],nn,K,byrow = F)*bsp)
    }
    
    # theta estimation
    
    invzz = try(MASS::ginv(t(ZZ[-c(1:mm),])%*%(ZZ[-c(1:mm),])),silent = T)
    
    if(class(invzz)[1] == "try-error"){
      invzz = g.inverse(t(ZZ[-c(1:mm),])%*%(ZZ[-c(1:mm),]))
    }
    
    theta_est = t(yy[-c(1:mm),])%*%ZZ[-c(1:mm),]%*%invzz
    
    loss = sum((yy[-c(1:mm),] - ZZ[-c(1:mm),]%*%t(theta_est))^2)
    
    return(loss/pp)
    
  }
  
  
  LnAlpha = function(alpha,theta,vv,yy,xx,K,mm,...){
    
    alpha = rescale(alpha)
    
    pp = NCOL(yy)
    dd = NCOL(vv)
    qq = NCOL(xx)
    nn = NROW(yy)
    
    ss = vv%*%alpha
    
    bsp = bs_new(ss,K,mm)
    
    ZZ = bsp
    
    for (ii in 1:qq) {
      ZZ = cbind(ZZ,matrix(xx[,ii],nn,K,byrow = F)*bsp)
    }
    
    loss = sum((yy[-c(1:mm),] - ZZ[-c(1:mm),]%*%t(theta))^2)
    
    return(loss/pp)
    
  }
  
  nn = NROW(yy)
  pp = NCOL(yy)
  
  if(is.null(zz)){
  zz = xx  
  }
  
  yy = as.matrix(yy)
  xx = as.matrix(xx)
  zz = as.matrix(zz)
  
  dd = NCOL(zz)
  qq = NCOL(xx)
    
  vv = matrix(NA,nrow = nn,ncol= dd)
  
  if(dd == 1){
    for(i in (mm+1):nn){if(mm==1){vv[i,] = zz[(i-1),]}else{vv[i,] = sum(zz[(i-1):(i-mm),])/mm}}
  }else{
    for(i in (mm+1):nn){if(mm==1){vv[i,] = zz[(i-1),]}else{vv[i,] = colSums(zz[(i-1):(i-mm),])/mm}}
    }
  
  #initialized 
  if(is.null(alpha_inl)){
    opt_inl = Rsolnp::gosolnp(fun = function(z){LnAlpha_inl(z,vv,yy,xx,K,mm)},LB = c(0,rep(LB,dd-1)),UB = rep(UB,dd), n.sim = n.sim, n.restarts = n.restarts,control = list(trace = 0))
    alpha_inl = rescale(opt_inl$pars) 
    obj_inl =  tail(opt_inl$values,1) 
  }else{
    
    if(length(alpha_inl) != dd){
      stop("The dimension of alpha doesn't match Z!")
    }
    
    obj_inl = 99999999
  }
  
  if(is.null(K)){
    
    bic = vector()
    alpha_tol = list()
    theta_tol = list()
    
    for (K in Kmin:Kmax){
      
      alpha = rep(0,dd) 
      alpha_new =  alpha_inl
      iter = 0
      eps = eps
      obj = vector()
      obj_t = obj_inl+1
      obj_t1 = obj_inl
      
      while( sum(abs(alpha - alpha_new)) > eps   ){  ##&& obj_t1 <  obj_t
        obj_t = obj_t1
        alpha = alpha_new
        
        ss = vv%*%alpha
        
        bsp = bs_new(ss,K,mm)
        
        theta = rep(1,(qq+1)*K)
        
        ZZ = bsp
        
        for (ii in 1:qq) {
          ZZ = cbind(ZZ,matrix(xx[,ii],nn,K,byrow = F)*bsp)
          
        }
        
        # theta estimation
        invzz = try(MASS::ginv(t(ZZ[-c(1:mm),])%*%(ZZ[-c(1:mm),])),silent = T)
        
        if(class(invzz)[1] == "try-error"){
          invzz = g.inverse(t(ZZ[-c(1:mm),])%*%(ZZ[-c(1:mm),]))
        }
        
        theta_est = t(yy[-c(1:mm),])%*%ZZ[-c(1:mm),]%*%invzz
        
        #alpha estimation
        
        opt = Rsolnp::solnp(alpha,fun = function(z){LnAlpha(z,theta_est,vv,yy,xx,K,mm)},control = list(trace = 0))
        
        alpha_new = rescale(opt$pars)
        obj_t1 = tail(opt$values,1)
        
        if(alpha_new[1] < 0){
          alpha_new = - alpha_new
        }
        
        iter = iter + 1
        
        obj[iter] <-  tail(opt$values,1)
        
        if(trace == 1){
        cat(paste0("########  ","K= ",K,"  ####  Obj:",sprintf("%0.3f",obj[iter]),"  ####  ","Alpha Diff:",sprintf("%0.4f",sum(abs(alpha - alpha_new))),"  ####  Iteration:", iter,"  ########"),"\n") 
        }
                
        
        if(iter > itermax-1){
          break
        }
      }
      
      theta_tol[[K-Kmin+1]] <- theta_est
      
      alpha_tol[[K-Kmin+1]] <- alpha_new
      
      bic[K-Kmin+1] = log(opt$value) + 0.5*K*log(nn)/nn 
      
      alpha_inl = alpha_new
    }
    
    if(Maxisbest == T){
      K_best = Kmax - Kmin + 1
    }else{
      K_best = which(bic == min(bic))
    }
    
    alpha_hat = alpha_tol[[K_best]]
    theta_hat = theta_tol[[K_best]]
    
    K_best = K_best + Kmin - 1
    
    
    
  }else{
    
    alpha_new = alpha_inl 
    alpha = rep(0,dd) 
    iter = 0
    eps = eps
    obj = vector()
    obj_t = obj_inl+1
    obj_t1 = obj_inl
    
    while( sum(abs(alpha - alpha_new)) > eps  ){ #&& obj_t1 <  obj_t
      
      obj_t = obj_t1
      alpha = alpha_new
      
      ss = vv%*%alpha
      
      bsp = bs_new(ss,K,mm)
      
      theta = rep(1,(qq+1)*K)
      
      ZZ = bsp
      
      for (ii in 1:qq) {
        ZZ = cbind(ZZ,matrix(xx[,ii],nn,K,byrow = F)*bsp)
        
      }
      
      # theta estimation
      
      invzz = try(MASS::ginv(t(ZZ[-c(1:mm),])%*%(ZZ[-c(1:mm),])),silent = T)
      
      if(class(invzz)[1] == "try-error"){
        invzz = g.inverse(t(ZZ[-c(1:mm),])%*%(ZZ[-c(1:mm),]))
      }
      
      theta_est = t(yy[-c(1:mm),])%*%ZZ[-c(1:mm),]%*%invzz
      
      #alpha estimation
      
      opt = Rsolnp::solnp(alpha,fun = function(z){LnAlpha(z,theta_est,vv,yy,xx,K,mm)},control = list(trace = 0))
      
      alpha_new = rescale(opt$pars)
      obj_t1 = tail(opt$values,1)
      
      if(alpha_new[1] < 0){
        alpha_new = - alpha_new
      }
      
      iter = iter + 1
      
      obj[iter] <-  tail(opt$values,1)
      
      if(trace == 1){
        cat(paste0("########  ","Obj:",sprintf("%0.3f",obj[iter]),"  ####  ","Alpha Diff:",sprintf("%0.4f",sum(abs(alpha - alpha_new))),"  ####  Iteration:", iter,"  ########"),"\n") 
      }
            
      if(iter > itermax-1){
        break
      }
      
    }
    
    K_best = K
    bic    = NULL
    alpha_tol = NULL
    theta_tol = NULL
    alpha_hat = alpha_new
    theta_hat = theta_est
  }
  
  #### fit and forecast
  
  ss_fit = vv%*%alpha_hat
  K = K_best
  
  if(mm > 1){
    vvv = matrix(NA,nrow = nn,ncol = dd)
    for(i in mm:nn){
      vvv[i,] = colSums(zz[(i):(i-mm+1),])/mm
    }
  }else{
    vvv = zz
  }
  
  ss_fore = vvv%*%alpha_hat
  
  ss_now = ss_fore[nn]
  
  
  ss_fit = na.omit(ss_fit) ######### na omit !!!!!!!!!!!!!!
 
  ss_down = quantile(ss_fit,trim) 
  ss_up   = quantile(ss_fit,1-trim)
  
  ss_range = which(ss_fit > ss_down & ss_fit < ss_up )
  
  bsp = bs_new(ss_fit,K,mm)
  
  mu_range = bsp[ss_range,]%*%t(theta_hat[,1:K])

  mu_fit = matrix(NA,nn,pp)
  mu_fore = vector()
  
  for (jjj in 1:pp) {
    mu_fit[-c(1:mm),jjj] = approx(ss_fit[ss_range] , mu_range[,jjj], ss_fit , rule = 2)$y 
    
    mu_fore[jjj] = approx(ss_fit[ss_range] , mu_range[,jjj], ss_now , rule = 2)$y 
  }
  
  
  ##### beta fit
  
  beta_fit = array(NA,dim = c(nn,pp,qq))
  
  beta_fore = matrix(NA,pp,qq)
  
  for (jj in 1:qq) {
    beta_range = bsp[ss_range,]%*%t(theta_hat[,(jj*K+1):((jj+1)*K)])
    
    for (ii in 1:pp){
      beta_fit[-c(1:mm),ii,jj] = approx(ss_fit[ss_range] , beta_range[,ii], ss_fit , rule = 2)$y
      
      beta_fore[ii,jj] = approx(ss_fit[ss_range], beta_range[,ii], ss_now , rule = 2)$y
      
    }
    
  }
  
  # error estimate
  
  
  eps_est  = matrix(NA,nn,pp)
  
  for (uuu in 1:pp) {
    
    eps_est[,uuu] = yy[,uuu] - mu_fit[,uuu]  -  rowSums(beta_fit[,uuu,]*xx) 
    
  }          
  
  con = list(alpha_est   =   alpha_hat,
             mu_est      =   mu_fore,
             beta_est    =   beta_fore,
             eps_est     =   eps_est,
             beta_fit    =   beta_fit,
             mu_fit      =   mu_fit,
             theta_est   =   theta_hat,
             K_best      =   K_best,
             alpha_tol   =   alpha_tol,
             theta_tol   =   theta_tol,
             bic         =   bic,
             vv          =   vv ,
             ss_now      =   ss_now,
             ss_fit      =   ss_fit
  )
  
  return(con)
  
  
}



beta.bspline = function(yy,xx,zz = NULL,alpha,mm = 1,K = 10,Kmin = 5, Kmax = 15, trim = 0.05,Maxisbest = F,...){
  
  nn = NROW(yy)
  qq = NCOL(xx)
  pp = NCOL(yy)
  
  if(is.null(zz)){
    zz = xx
  }
  
  yy = as.matrix(yy)
  xx = as.matrix(xx)
  zz = as.matrix(zz)
  
  dd = NCOL(zz)
  
  vv = matrix(NA,nrow=nn,ncol= dd)
  
  if(dd == 1){
    for(i in (mm+1):nn){if(mm==1){vv[i,] = xx[(i-1),]}else{vv[i,] = sum(xx[(i-1):(i-mm),])/mm}}
  }else{
    for(i in (mm+1):nn){if(mm==1){vv[i,] = xx[(i-1),]}else{vv[i,] = colSums(xx[(i-1):(i-mm),])/mm}}}
  
  
  if(is.null(K)){
    
    bic = vector()
    theta_tol = list()
    for (K in Kmin:Kmax){
      
      ss = vv%*%alpha
      
      bsp = bs_new(ss,K,mm)
      
      theta = rep(1,(qq+1)*K)
      
      ZZ = bsp
      
      for (ii in 1:qq) {
        ZZ = cbind(ZZ,matrix(xx[,ii],nn,K,byrow = F)*bsp)
        
      }
      
      # theta estimation
      
      invzz = try(MASS::ginv(t(ZZ[-c(1:mm),])%*%(ZZ[-c(1:mm),])),silent = T)
      
      if(class(invzz)[1] == "try-error"){
        invzz = g.inverse(t(ZZ[-c(1:mm),])%*%(ZZ[-c(1:mm),]))
      }
      
      theta_est = t(yy[-c(1:mm),])%*%ZZ[-c(1:mm),]%*%invzz
      
      obj <- sum((yy[-c(1:mm),] - ZZ[-c(1:mm),]%*%t(theta_est))^2)/pp
      
      theta_tol[[K-Kmin+1]] <- theta_est
      
      bic[K-Kmin+1] = log(obj) + 0.5*K*log(nn)/nn 
      
    }
    
    if(Maxisbest == T){
      K_best = Kmax - Kmin + 1
    }else{
      K_best = which(bic == min(bic))
    }
    
    theta_hat = theta_tol[[K_best]]
    
    K_best = K_best + Kmin - 1
    
  }else{
    
    ss = vv%*%alpha
    
    bsp = bs_new(ss,K,mm)
    
    theta = rep(1,(qq+1)*K)
    
    ZZ = bsp
    
    for (ii in 1:qq) {
      ZZ = cbind(ZZ,matrix(xx[,ii],nn,K,byrow = F)*bsp)
      
    }
    
    # theta estimation
    
    invzz = try(MASS::ginv(t(ZZ[-c(1:mm),])%*%(ZZ[-c(1:mm),])),silent = T)
    
    if(class(invzz)[1] == "try-error"){
      invzz = g.inverse(t(ZZ[-c(1:mm),])%*%(ZZ[-c(1:mm),]))
    }
    
    theta_est = t(yy[-c(1:mm),])%*%ZZ[-c(1:mm),]%*%invzz
    
    
    K_best = K
    bic    = NULL
    theta_tol = NULL
    theta_hat = theta_est
  }
  
  #### fit and forecast
  
  ss_fit = vv%*%alpha
  K = K_best
  
  if(mm > 1){
    vvv = matrix(NA,nrow = nn,ncol = dd)
    for(i in mm:nn){
      vvv[i,] = colSums(zz[(i):(i-mm+1),])/mm
    }
  }else{
    vvv = zz
  }
  
  ss_fore = vvv%*%alpha
  
  ss_now = ss_fore[nn]
  
  
  ss_fit = na.omit(ss_fit) ######### na omit !!!!!!!!!!!!!!

  ss_down = quantile(ss_fit,trim) 
  ss_up   = quantile(ss_fit,1-trim)
  
  ss_range = which(ss_fit > ss_down & ss_fit < ss_up )
  
  bsp = bs_new(ss_fit,K,mm)
  
  mu_range = bsp[ss_range,]%*%t(theta_hat[,1:K])
  
  mu_fit = matrix(NA,nn,pp)
  mu_fore = vector()
  
  for (jjj in 1:pp) {
    mu_fit[-c(1:mm),jjj] = approx(ss_fit[ss_range] , mu_range[,jjj], ss_fit , rule = 2)$y 
    
    mu_fore[jjj] = approx(ss_fit[ss_range] , mu_range[,jjj], ss_now , rule = 2)$y 
  }
  
  
  ##### beta fit
  
  beta_fit = array(NA,dim = c(nn,pp,qq))
  
  beta_fore = matrix(NA,pp,qq)
  
  for (jj in 1:qq) {
    beta_range = bsp[ss_range,]%*%t(theta_hat[,(jj*K+1):((jj+1)*K)])
    
    for (ii in 1:pp){
      beta_fit[-c(1:mm),ii,jj] = approx(ss_fit[ss_range] , beta_range[,ii], ss_fit , rule = 2)$y
      
      beta_fore[ii,jj] = approx(ss_fit[ss_range], beta_range[,ii], ss_now , rule = 2)$y
      
    }
    
  }
  
  # error estimate
  
  eps_est  = matrix(NA,nn,pp)
  
  for (uuu in 1:pp) {
    
    eps_est[,uuu] = yy[,uuu] - mu_fit[,uuu]  -  rowSums(beta_fit[,uuu,]*xx) 
    
  }          
  
  con = list(alpha_est   =   alpha,
             mu_est      =   mu_fore,
             beta_est    =   beta_fore,
             eps_est     =   eps_est,
             beta_fit    =   beta_fit,
             mu_fit      =   mu_fit,
             theta_est   =   theta_hat,
             K_best      =   K_best,
             theta_tol   =   theta_tol,
             bic         =   bic,
             vv          =   vv ,
             ss_now      =   ss_now,
             ss_fit      =   ss_fit
  )
  
  return(con)
  
  
}



bs <- splines::bs

GLR.test = function(yy,xx,mm = 1,KK = 10,rep = 100,...){
  
  yy = as.matrix(yy)
  xx = as.matrix(xx)
  
  n = NCOL(yy)
  t = NROW(yy)
  
  lmInitial = lm(yy ~ xx)
  
  #(1)BSPLINE-ECSNP estimation
  alphaInitial = rescale(rowMeans(lmInitial$coefficients)[-1])
  
  # alphaInitial_fin = alpha_est_inl(yy,xx,mm = mm,alpha_inl = alphaInitial,itermax = 50)
  
  result_bsp = coef.bspline(yy,xx,mm = mm,alpha_inl = alphaInitial,K = KK,trim = 0.01,itermax = 1000)
  
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
    
    result_bsp_bs = beta.bspline(yy_bs,xx,mm = mm,alpha = alpha_est,K = KK,trim = 0.01)
    
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


coef.bspline = function(yy,xx,alpha_inl,mm = 1,K = NULL, Kmax = 5, eps = 10^-5, itermax = 100, num_gp = 500, trim = 0.05,Maxisbest = F,...){
  
  LnAlpha = function(alpha,theta,vv,yy,xx,K,mm,...){
    
    alpha = rescale(alpha)
    
    ss = vv%*%alpha
    
    bsp = bs_new(ss,K,mm)
    
    ZZ = bsp
    
    for (ii in 1:dd) {
      ZZ = cbind(ZZ,matrix(xx[,ii],nn,K,byrow = F)*bsp)
    }
    
    loss = sum((yy[-c(1:mm),] - ZZ[-c(1:mm),]%*%t(theta))^2)
    
    return(loss)
    
  }
  
  nn = NROW(yy)
  dd = NCOL(xx)
  pp = NCOL(yy)
  
  vv = matrix(NA,nrow=nn,ncol= dd)
  
  if(dd == 1){
    for(i in (mm+1):nn){if(mm==1){vv[i,] = xx[(i-1),]}else{vv[i,] = sum(xx[(i-1):(i-mm),])/mm}}
  }else{
    for(i in (mm+1):nn){if(mm==1){vv[i,] = xx[(i-1),]}else{vv[i,] = colSums(xx[(i-1):(i-mm),])/mm}}}
  
  
  if(is.null(K)){
    
    bic = vector()
    alpha_tol = list()
    theta_tol = list()
    for (K in 3:Kmax){
      
      alpha = rep(0,dd) 
      alpha_new =  alpha_inl
      iter = 0
      eps = eps
      obj = vector()
      
      while( sum(abs(alpha - alpha_new)) > eps  ){
        
        alpha = alpha_new
        
        ss = vv%*%alpha
        
        bsp = bs_new(ss,K,mm)
        
        theta = rep(1,(dd+1)*K)
        
        ZZ = bsp
        
        for (ii in 1:dd) {
          ZZ = cbind(ZZ,matrix(xx[,ii],nn,K,byrow = F)*bsp)
          
        }
        
        # theta estimation
        
        theta_est = t(yy[-c(1:mm),])%*%ZZ[-c(1:mm),]%*%MASS::ginv(t(ZZ[-c(1:mm),])%*%(ZZ[-c(1:mm),]))
        
        #alpha estimation
        
        opt = optim(par = alpha, fn = function(z){LnAlpha(z,theta_est,vv,yy,xx,K,mm)}, method = "L-BFGS-B")
        
        alpha_new = rescale(opt$par) 
        
        if(alpha_new[1] < 0){
          alpha_new = - alpha_new
        }
        
        iter = iter + 1
        
        obj[iter] <- opt$value/pp
        
        print(paste0("K=",K," ",obj[iter]))
        
        if(iter > itermax){
          break
        }
      }
      
      theta_tol[[K-2]] <- theta_est
      
      alpha_tol[[K-2]] <- alpha_new
      
      bic[K-2] = log(opt$value) + 0.5*K*log(nn)/nn 
      
      alpha_inl = alpha_new
    }
    
    if(Maxisbest == T){
      K_best = Kmax
    }else{
      K_best = which(bic == min(bic))
    }
    
    alpha_hat = alpha_tol[[K_best]]
    theta_hat = theta_tol[[K_best]]
    
    K_best = K_best + 2
    
    
    
  }else{
    
    alpha = rep(0,dd) 
    alpha_new =  alpha_inl
    iter = 0
    eps = eps
    obj = vector()
    
    while( sum(abs(alpha - alpha_new)) > eps  ){
      
      alpha = alpha_new
      
      ss = vv%*%alpha
      
      bsp = bs_new(ss,K,mm)
      
      theta = rep(1,(dd+1)*K)
      
      ZZ = bsp
      
      for (ii in 1:dd) {
        ZZ = cbind(ZZ,matrix(xx[,ii],nn,K,byrow = F)*bsp)
        
      }
      
      # theta estimation
      
      theta_est = t(yy[-c(1:mm),])%*%ZZ[-c(1:mm),]%*%MASS::ginv(t(ZZ[-c(1:mm),])%*%(ZZ[-c(1:mm),]))
      
      #alpha estimation
      
      opt = optim(par = alpha, fn = function(z){LnAlpha(z,theta_est,vv,yy,xx,K,mm)}, method = "L-BFGS-B")
      
      alpha_new = rescale(opt$par) 
      
      if(alpha_new[1] < 0){
        alpha_new = - alpha_new
      }
      
      iter = iter + 1
      
      obj[iter] <- opt$value/pp
      
      print(obj[iter])
      
      if(iter > itermax){
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
      vvv[i,] = colSums(xx[(i):(i-mm+1),])/mm
    }
  }else{
    vvv = xx
  }
  
  ss_fore = vvv%*%alpha_hat
  
  ss_now = ss_fore[nn]
  
  
  ss_fit = na.omit(ss_fit) ######### na omit !!!!!!!!!!!!!!

  # zz = calc_grid_points(na.omit(vv),alpha_hat,num_gp = num_gp) #(not good) first, give a grid net; second, delete tail value; third, use interpolate give the fit beta and mu
  
  ss_down = quantile(ss_fit,trim) 
  ss_up   = quantile(ss_fit,1-trim)
  
  ss_range = which(ss_fit > ss_down & ss_fit < ss_up )
  
  bsp = bs_new(ss_fit,K,mm)
  
  mu_range = bsp[ss_range,]%*%t(theta_hat[,1:K])
  

  # bsp = bs_new(zz,K,mm)
  # 
  # zz_range = seq(from = round(num_gp*trim),to = (num_gp*(1-trim)))
  # 
  # mu_range = bsp[zz_range,]%*%t(theta_est[,1:K])
  
  mu_fit = matrix(NA,nn,pp)
  mu_fore = vector()
    
  for (jjj in 1:pp) {
    mu_fit[-c(1:mm),jjj] = approx(ss_fit[ss_range] , mu_range[,jjj], ss_fit , rule = 2)$y 
    
    mu_fore[jjj] = approx(ss_fit[ss_range] , mu_range[,jjj], ss_now , rule = 2)$y 
  }
  
 
  ##### beta fit
  
  beta_fit = array(NA,dim = c(nn,pp,dd))
  
  beta_fore = matrix(NA,pp,dd)
  
  for (jj in 1:dd) {
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



beta.bspline = function(yy,xx,alpha,mm = 1,K = NULL, Kmax = 15, num_gp = 500, trim = 0.025,Maxisbest = F,...){

  nn = NROW(yy)
  dd = NCOL(xx)
  pp = NCOL(yy)
  
  vv = matrix(NA,nrow=nn,ncol= dd)
  
  if(dd == 1){
    for(i in (mm+1):nn){if(mm==1){vv[i,] = xx[(i-1),]}else{vv[i,] = sum(xx[(i-1):(i-mm),])/mm}}
  }else{
    for(i in (mm+1):nn){if(mm==1){vv[i,] = xx[(i-1),]}else{vv[i,] = colSums(xx[(i-1):(i-mm),])/mm}}}
  
  
  if(is.null(K)){
    
    bic = vector()
    theta_tol = list()
    for (K in 3:Kmax){
      
        ss = vv%*%alpha
        
        bsp = bs_new(ss,K,mm)
        
        theta = rep(1,(dd+1)*K)
        
        ZZ = bsp
        
        for (ii in 1:dd) {
          ZZ = cbind(ZZ,matrix(xx[,ii],nn,K,byrow = F)*bsp)
          
        }
        
        # theta estimation
        
      theta_est = t(yy[-c(1:mm),])%*%ZZ[-c(1:mm),]%*%MASS::ginv(t(ZZ[-c(1:mm),])%*%(ZZ[-c(1:mm),]))

      obj <- sum((yy[-c(1:mm),] - ZZ[-c(1:mm),]%*%t(theta_est))^2)/pp
        
      theta_tol[[K-2]] <- theta_est
  
      bic[K-2] = log(obj) + 0.5*K*log(nn)/nn 
      
    }
    
    if(Maxisbest == T){
      K_best = Kmax - 2
    }else{
      K_best = which(bic == min(bic))
    }

    theta_hat = theta_tol[[K_best]]
    
    K_best = K_best + 2

  }else{
    
      ss = vv%*%alpha
      
      bsp = bs_new(ss,K,mm)
      
      theta = rep(1,(dd+1)*K)
      
      ZZ = bsp
      
      for (ii in 1:dd) {
        ZZ = cbind(ZZ,matrix(xx[,ii],nn,K,byrow = F)*bsp)
        
      }
      
      # theta estimation
      
      theta_est = t(yy[-c(1:mm),])%*%ZZ[-c(1:mm),]%*%MASS::ginv(t(ZZ[-c(1:mm),])%*%(ZZ[-c(1:mm),]))
      

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
      vvv[i,] = colSums(xx[(i):(i-mm+1),])/mm
    }
  }else{
    vvv = xx
  }
  
  ss_fore = vvv%*%alpha
  
  ss_now = ss_fore[nn]
  
  
  ss_fit = na.omit(ss_fit) ######### na omit !!!!!!!!!!!!!!
  
  # zz = calc_grid_points(na.omit(vv),alpha_hat,num_gp = num_gp) #(not good) first, give a grid net; second, delete tail value; third, use interpolate give the fit beta and mu
  
  ss_down = quantile(ss_fit,trim) 
  ss_up   = quantile(ss_fit,1-trim)
  
  ss_range = which(ss_fit > ss_down & ss_fit < ss_up )
  
  bsp = bs_new(ss_fit,K,mm)
  
  mu_range = bsp[ss_range,]%*%t(theta_hat[,1:K])
  
  
  # bsp = bs_new(zz,K,mm)
  # 
  # zz_range = seq(from = round(num_gp*trim),to = (num_gp*(1-trim)))
  # 
  # mu_range = bsp[zz_range,]%*%t(theta_est[,1:K])
  
  mu_fit = matrix(NA,nn,pp)
  mu_fore = vector()
  
  for (jjj in 1:pp) {
    mu_fit[-c(1:mm),jjj] = approx(ss_fit[ss_range] , mu_range[,jjj], ss_fit , rule = 2)$y 
    
    mu_fore[jjj] = approx(ss_fit[ss_range] , mu_range[,jjj], ss_now , rule = 2)$y 
  }
  
  
  ##### beta fit
  
  beta_fit = array(NA,dim = c(nn,pp,dd))
  
  beta_fore = matrix(NA,pp,dd)
  
  for (jj in 1:dd) {
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



face_roll_bspline = function(x,KK,...){
  
  tt_year = x
  
  data_sample_i = data_sample[(tt_year):(tt_year+TT_train-1),]
  
  if(tt_year == tail(TT_test_bk,1)){
    return_sample_i = data_sample[(tt_year+TT_train):TT_tol,]
    name_i = row.names(data_sample)[(tt_year+TT_train):TT_tol]
  }else{
    return_sample_i = data_sample[(tt_year+TT_train):(tt_year+TT_train+bk-1),]
    name_i = row.names(data_sample)[(tt_year+TT_train):(tt_year+TT_train+bk-1)]
  }
  
  
  ff_sample_i = as.matrix(ff_sample[(tt_year):(tt_year+TT_train-1)])
  ff3_sample_i = as.matrix(ff3_sample[(tt_year):(tt_year+TT_train-1),])
  
  ############### 2.1 estimation procedure ################
  
  lmInitial = lm(data_sample_i ~ ff3_sample_i)
  
  #(1)FACE-ECSNP estimation
  alphaInitial = rescale(rowMeans(lmInitial$coefficients)[-1])
  
  alphaInitial_fin = alpha_est_inl(data_sample_i,ff3_sample_i,alpha_inl = alphaInitial)
  
  result_bsp = coef.bspline(data_sample_i,ff3_sample_i,alpha_inl = alphaInitial_fin$alpha,K = KK,trim = 0.025,itermax = 500)
  
    
  # result_bsp$K_best
  
  # plot(result_bsp$bic)
   
  alpha_est =  result_bsp$alpha_est
  u_est     =  result_bsp$mu_est
  beta_est  =  result_bsp$beta_est
  K_best    =  result_bsp$K_best
  
  beta_fit = result_bsp$beta_fit
  
  eps_est  = result_bsp$eps_est
  
  
  result_face = result_roll_face[[tt_year]]
  
  factor_moments  =  result_roll_face[[tt_year]]$factor_moments
  
  factor_moments_ica = factor_moments$ica_moments
  factor_moments_dcc = factor_moments$dcc_moments
  factor_moments_copula = factor_moments$copula_moments
  
  vff <- vars::VAR(ff3_sample_i,type = "const")
  f_fore = predict(vff,n.ahead = 1)
  mu_f3_fore = c(f_fore$fcst[[1]][1],f_fore$fcst[[2]][1],f_fore$fcst[[3]][1])
  
  #########eps moments
  
  e_var_fore_mf_con <- vector()
  e_skew_fore_mf_con <- vector()
  e_kurt_fore_mf_con <- vector()
  
  con_residuals <- list()
  
  for (kk in 1:n){
    e_ik = eps_est[-(1:mm),kk]
    
    est_ei_snp =  TGC_est(e_ik,var.model = "sGARCH", var.targeting = F, 
                                    var.distribution = "sged", tgc.type = "leverage", 
                                    tgc.targeting = F, mean.model = list(armaOrder = c(0,0)), CTGC = T, rep_sim = 10)
    
    con_residuals[[kk]] <- est_ei_snp
    
    e_var_fore_mf_con[kk] <- est_ei_snp$result_moment$mm.fore[2]
    e_skew_fore_mf_con[kk] <- est_ei_snp$result_moment$mm.fore[3]
    e_kurt_fore_mf_con[kk] <- est_ei_snp$result_moment$mm.fore[4]
    
  }
  
  face_residual_moments = list(e_var_fore_mf_con,e_skew_fore_mf_con,e_kurt_fore_mf_con) 
  
  #save result of vc-mf-tvsnp
  result_roll_bs            <- list(alpha            =  alpha_est, 
                                      mu             =  u_est, 
                                      mu_f           =  mu_f3_fore,
                                      beta           =  beta_est, 
                                      factor_moments =  factor_moments, 
                                      eps_moments    =  face_residual_moments)
  
  ###portfolio
  #2.1 MZ-3F
  beta_3f = t(lmInitial$coefficients[-1,])
  mu_3f = lmInitial$coefficients[1,] + beta_3f%*%colMeans(ff3_sample_i)
  
  #3.1 SEMI-3F-ECSNP-ICA
  beta_3f_face  = beta_est
  mu_3f_face = u_est + beta_3f_face%*%mu_f3_fore
  # mu_3f_face = colMeans(data_sample_i)
  mmf_3f_snp_ica = factor_moments_ica
  mme_3f_snp_face = face_residual_moments
  
  EU_FACEECSNP_3F_ica = function(z){
    mmP = Portfolio.Cumulants.Mat(z,mm_factor = mmf_3f_snp_ica,mm_eps = mme_3f_snp_face,A = beta_3f_face)
    obj = Obj.EU(mmP,gamma = gamma);return(obj)
  }
  eq_3f_face = function(x){z1=sum(x);return(z1)}
  ineq_3f_face = function(x){sum(x*mu_3f_face)}
  ineq_3f = function(x){sum(x*mu_3f)}
  
  #3.2 SEMI-3F-ECSNP-DCC
  mmf_3f_snp_dcc = factor_moments_dcc
  
  EU_FACEECSNP_3F_dcc = function(z){
    mmP = Portfolio.Cumulants.Mat(z,mm_factor = mmf_3f_snp_dcc,mm_eps = mme_3f_snp_face,A = beta_3f_face)
    obj = Obj.EU(mmP,gamma = gamma);return(obj)
  }
  
  #3.3 SEMI-3F-ECSNP-copula
  mmf_3f_snp_copula = factor_moments_copula
  
  EU_FACEECSNP_3F_copula = function(z){
    mmP = Portfolio.Cumulants.Mat(z,mm_factor = mmf_3f_snp_copula,mm_eps = mme_3f_snp_face,A = beta_3f_face)
    obj = Obj.EU(mmP,gamma = gamma);return(obj)
  }
  
  
  lb = -5;ub = 5;delta = 1;gamma = 5
  ############################2.3 portfolio selection###########################
 
  rsol_3f_snp_ica    = Rsolnp::solnp(pars = rep(1/n,n),fun = EU_FACEECSNP_3F_ica,eqfun = eq_3f_face,ineqfun = ineq_3f,ineqLB = delta,ineqUB = 100,eqB = 1,LB = rep(lb,n),UB = rep(ub,n))
  rsol_3f_snp_dcc    = Rsolnp::solnp(pars = rep(1/n,n),fun = EU_FACEECSNP_3F_dcc,eqfun = eq_3f_face,ineqfun = ineq_3f,ineqLB = delta,ineqUB = 100,eqB = 1,LB = rep(lb,n),UB = rep(ub,n))
  rsol_3f_snp_copula = Rsolnp::solnp(pars = rep(1/n,n),fun = EU_FACEECSNP_3F_copula,eqfun = eq_3f_face,ineqfun = ineq_3f,ineqLB = delta,ineqUB = 100,eqB = 1,LB = rep(lb,n),UB = rep(ub,n))
  
  rsol_3f_face_ica    = Rsolnp::solnp(pars = rep(1/n,n),fun = EU_FACEECSNP_3F_ica,eqfun = eq_3f_face,ineqfun = ineq_3f_face,ineqLB = delta,ineqUB = 100,eqB = 1,LB = rep(lb,n),UB = rep(ub,n))
  rsol_3f_face_dcc    = Rsolnp::solnp(pars = rep(1/n,n),fun = EU_FACEECSNP_3F_dcc,eqfun = eq_3f_face,ineqfun = ineq_3f_face,ineqLB = delta,ineqUB = 100,eqB = 1,LB = rep(lb,n),UB = rep(ub,n))
  rsol_3f_face_copula = Rsolnp::solnp(pars = rep(1/n,n),fun = EU_FACEECSNP_3F_copula,eqfun = eq_3f_face,ineqfun = ineq_3f_face,ineqLB = delta,ineqUB = 100,eqB = 1,LB = rep(lb,n),UB = rep(ub,n))

  
  
  w_3f_snp_ica       =    rsol_3f_snp_ica$pars
  w_3f_snp_dcc       =    rsol_3f_snp_dcc$pars
  w_3f_snp_copula    =    rsol_3f_snp_copula$pars
  
  
  w_3f_face_ica       =    rsol_3f_face_ica$pars
  w_3f_face_dcc       =    rsol_3f_face_dcc$pars
  w_3f_face_copula    =    rsol_3f_face_copula$pars
  
  
  eq = rep(1/n,n)
  
  u_total = cbind(eq, 
                  f3_snp_ica        =  w_3f_snp_ica,
                  f3_snp_dcc        =  w_3f_snp_dcc,
                  f3_snp_copula     =  w_3f_snp_copula,
                  f3_face_ica       =  w_3f_face_ica,
                  f3_face_dcc       =  w_3f_face_dcc,
                  f3_face_copula    =  w_3f_face_copula
  )
  
if(NROW(return_sample_i) == 1){
    r_total_u =  t(return_sample_i)%*%u_total 
  }else{
    r_total_u =  return_sample_i%*%u_total
}
  
  con = list(r_total_u = r_total_u,u_total = u_total,K_best = K_best,beta_est = beta_est , alpha_est = alpha_est, result_roll_bs = result_roll_bs)
  
  return(con)
}

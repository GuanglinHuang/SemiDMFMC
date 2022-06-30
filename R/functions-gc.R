#' @useDynLib SmeiDMFMC
#' @importFrom Rcpp sourceCpp

ugarchforecast <- rugarch::ugarchforecast

GC_est = function(ff,var.model = 'sGARCH',var.targeting = F,var.distribution = 'sstd',
                  gc.type  = "leverage",gc.targeting = F,mean.model = list(armaOrder = c(1, 1)),CGC = FALSE){
  
  type  = gc.type
  
  spec = rugarch::ugarchspec(mean.model = mean.model,variance.model = list(model = var.model, garchOrder = c(1, 1),variance.targeting = var.targeting), distribution = var.distribution)
  
  fit_variance = rugarch::ugarchfit(spec, ff, solver = 'hybrid')
  
  z = fit_variance@fit$z  #residuals is ht*et, z is et,fitted.value is conditional mean
  
  hz = fit_variance@fit$residuals
  
  var_x_fit = fit_variance@fit$var
  var_x_s = hz^2
  mu_x_fit = fit_variance@fit$fitted.values
  
  coef_mv = fit_variance@fit$coef
  
  if(CGC == T){
    #$constant theta estimation$
    theta_con_inl = c(0,5)
    
    constant_est = optim(par = theta_con_inl, function(theta){
      con = likelihood_gc_con(z,theta)
      return(-con[[1]])
    },hessian = T)
    
    theta_est = constant_est$par
    std = sqrt(diag(solve(constant_est$hessian)))
    t_stat = constant_est$par/std
    p_value = 2*(1 - pnorm(abs(t_stat)))
    
    constant_con = cbind(theta_est,std,t_stat,p_value) #constant parameters estimation #return
    row.names(constant_con) <- c("theta1","theta2")
    theta_est = constant_con[,1] #return
    
    #moment forecasting
    theta_tv_inl = c(theta_est[1],rep(0,4),theta_est[2],rep(0,4))
    
    result_moment = gc_predict(fit_variance,ff,theta_tv_inl,type = type)
    
    result_con = list(
      snp.cof = theta_est, constant_con = constant_con,moment_fore = result_moment
    )
    
    return(result_con)
  }else{
    #$constant theta estimation$
    theta_con_inl = c(0,5)
    
    constant_est = optim(par = theta_con_inl, function(theta){
      con = likelihood_gc_con(z,theta)
      return(-con[[1]])
    },hessian = T)
    
    theta_est = constant_est$par
    std = sqrt(diag(solve(constant_est$hessian)))
    t_stat = constant_est$par/std
    p_value = 2*(1 - pnorm(abs(t_stat)))
    
    constant_con = cbind(theta_est,std,t_stat,p_value) #constant parameters estimation #return
    row.names(constant_con) <- c("theta1","theta2")
    theta_est = constant_con[,1] #return
    
    #moment forecasting - CSNP
    theta_tv_inl = c(theta_est[1],rep(0,4),theta_est[2],rep(0,4))
    
    result_moment = gc_predict(fit_variance,ff,theta_tv_inl,type = type)
    
    result_con = list(
      snp.cof = theta_est, constant_con = constant_con,moment_fore = result_moment
    )
    
    #$time varying theta estimation$
    if(!gc.targeting){
      ll_tol = vector()
      con_tol = list()
      theta_est_tol = list()
      for (bn in 1:rep_sim) {
        theta_tv_inl = c(0,runif(1,-0.8,0.8),runif(3,-0.3,0.3),5,runif(1,-0.8,0.8),runif(3,-0.3,0.3))
        con_tv = optim(par = theta_tv_inl,
                       function(s){theta_tv = s;
                       con = likelihood_gc_tv_rcpp(z,theta_tv,type = type);
                       return(-con)}, method = "BFGS", hessian = T)
        theta_tv_est = con_tv$par #return
        ll_tol[bn] <- likelihood_gc_tv_rcpp(z,theta_tv_est,type = type)
        con_tol[[bn]] <- con_tv
        theta_est_tol[[bn]] = theta_tv_est
      }
      
      pl = which(ll_tol == max(ll_tol))[1]
      con_tv = con_tol[[pl]]
      theta_tv_est = theta_est_tol[[pl]]
      
      std_tv = solve_hess(con_tv$hessian)
      t_stat_tv = con_tv$par/std_tv
      p_value_tv = 2*(1-pnorm(abs(t_stat_tv)))
      tv_con = cbind(theta_tv_est,std_tv,t_stat_tv,p_value_tv)#tv parameters estimation
      row.names(tv_con)<-c("omega1","alpha1","beta11","beta12","gamma1","omega2","alpha2","beta21","beta22","gamma2")
    }
    # targeting estimation
    if(gc.targeting){
      theta_tv_inl_tar = rep(0,8) #not finished
      con_tv_tar = optim(par = theta_tv_inl_tar, function(theta_tv){con = likelihood_gc_tv_targeting_rcpp(z,theta_tv,theta0 = theta_est,type = type);return(-con)}, method = "BFGS", hessian = T)
      theta_tv_est_tar = con_tv_tar$par
      std_tv_tar = solve_hess(con_tv_tar$hessian)
      t_stat_tv_tar = con_tv_tar$par/std_tv_tar
      p_value_tv_tar = 2*(1-pnorm(abs(t_stat_tv_tar)))
      tv_con_tar = cbind(theta_tv_est_tar,std_tv_tar,t_stat_tv_tar,p_value_tv_tar)#tv parameters estimation
      tv_con = rbind(constant_con[1,],tv_con_tar[1:4,],constant_con[2,],tv_con_tar[5:8,]) #return
      row.names(tv_con_tar)<-c("theta1","alpha1","beta11","beta12","gamma1","theta2","alpha2","beta21","beta22","gamma2")
      
      theta_tv_est = c(theta_est[1],theta_tv_est_tar[1:4],theta_est[2],theta_tv_est_tar[5:8]) #return
    }
    
    
    result_tv = list(
      gc.cof = theta_tv_est, tv_con = tv_con
    )
    
    #$moment forecasting$
    mm_fit = gc_fit_moment(z,theta_tv_est) #return
    
    skew_fit = mm_fit[,3]
    kurt_fit = mm_fit[,4]
    
    mm_fore = tgc_predict(fit_variance,ff,theta_tv_est,type = type)
    
    result_moment = list(
      mm.fit = mm_fit,
      mm.fore = mm_fore
    )
    return(list(result_con = result_con,
                result_tv = result_tv,
                result_moment = result_moment,
                fit_variance =fit_variance)
    )
  }
  
}


likelihood_gc_tv = function(z,theta_tv,type = c("n-leverage","leverage","linear"),..){
  t = length(z)
  if(type == "linear"){  #time-varying likelihood
    theta10_tv = theta_tv[1]
    theta11_tv = theta_tv[2]
    theta121_tv = theta_tv[3]
    theta122_tv = theta121_tv
    theta13_tv = 0

    theta20_tv = theta_tv[6]
    theta21_tv = theta_tv[7]
    theta221_tv = theta_tv[8]
    theta222_tv = theta221_tv
    theta23_tv = 0
  }

  if(type == "leverage"){
    theta10_tv = theta_tv[1]
    theta11_tv = theta_tv[2]
    theta121_tv = theta_tv[3]
    theta122_tv = theta_tv[4]
    theta13_tv = 0

    theta20_tv = theta_tv[6]
    theta21_tv = theta_tv[7]
    theta221_tv = theta_tv[8]
    theta222_tv = theta_tv[9]
    theta23_tv = 0
  }

  if(type == "n-leverage"){
    theta10_tv = theta_tv[1]
    theta11_tv = theta_tv[2]
    theta121_tv = theta_tv[3]
    theta122_tv = theta_tv[4]
    theta13_tv = theta_tv[5]

    theta20_tv = theta_tv[6]
    theta21_tv = theta_tv[7]
    theta221_tv = theta_tv[8]
    theta222_tv = theta_tv[9]
    theta23_tv = theta_tv[10]
  }

  theta1_t1 = 0
  theta2_t1 = 0

  theta1_tv = vector()
  theta2_tv = vector()

  #theta10 = theta0[1]
  #theta20 = theta0[2]

  theta1_tv[1] = theta10_tv
  theta2_tv[1] = theta20_tv


  lt_tol = vector()

  for (tt in 2:t) {
    zt1 = z[tt-1]
    zt = z[tt]

    theta1_t = theta10_tv + theta11_tv*theta1_t1 + theta121_tv*(1+ theta13_tv*abs(zt1))*max(zt1^3,0) + theta122_tv*(1+ theta13_tv*abs(zt1))*min(zt1^3,0)
    theta2_t = theta20_tv + theta21_tv*theta2_t1 + theta221_tv*(1+ theta23_tv*abs(zt1))*max(zt1^4,0) + theta222_tv*(1+ theta23_tv*abs(zt1))*min(zt1^4,0)

    gam1 = theta1_t/sqrt(6)
    gam2 = (theta2_t-3)/sqrt(24)

    lam = 1/(1+gam1^2+gam2^2)

    x_t = zt

    lt = log(lam) - (1/2)*x_t^2 +  2*log(abs(phi(x_t,theta1_t,theta2_t)))

    lt_tol = c(lt_tol,lt)

    theta1_tv[tt] = theta1_t
    theta2_tv[tt] = theta2_t

    theta1_t1 = theta1_t
    theta2_t1 = theta2_t
  }

  LL = sum(lt_tol)

  con = list(Lnf = LL,Theta = cbind(theta1_tv,theta2_tv))

  return(con)
}

likelihood_gc_con = function(z,theta,..){  #constant likehihood

  t = length(z)

  theta1_t = theta[1]
  theta2_t = theta[2]

  gam1 = theta1_t^2/6
  gam2 = (theta2_t-3)^2/24

  lam = 1/(1+gam1+gam2)

  x_t = z

  lt = log(lam) - (1/2)*x_t^2 +  2*log(abs(phi_gc(x_t,theta1_t,theta2_t)))

  LL = sum(lt)

  con = list(Lnf = LL,Theta = c(theta1_t,theta2_t),lt)

  return(con)
}

likelihood_gc_tv_rcpp = function(z,theta_tv,type = c("n-leverage","leverage","linear"),..){
  t = length(z)
  if(type == "linear"){  #time-varying likelihood
    theta_tv[5] <- 0
    theta_tv[4] <- theta_tv[3]
    theta_tv[9] <- theta_tv[8]
    theta_tv[10] <- 0
  }

  if(type == "leverage"){
    theta_tv[5] <- 0
    theta_tv[10] <- 0
  }

  if(type == "n-leverage"){
    theta_tv = theta_tv
  }
  con =lnf_gc_tv(z,theta_tv,length(z))

  return(con)
}

phi_gc = function(x,theta1,theta2,...){
  H3 = (x^3-3*x)/sqrt(6)
  H4 = (x^4-6*x^2+3)/sqrt(24)

  gam1 = theta1/sqrt(6)
  gam2 = (theta2-3)/sqrt(24)

  phi = 1 + gam1*H3 + gam2*H4
}

mm_gc_z = function(theta){
  theta1 = theta[,1]
  theta2 = theta[,2]

  m1_fit = 0
  m2_fit = 1
  m3_fit = theta1
  m4_fit = theta2-3

  mm = cbind(m1_fit,m2_fit,m3_fit,m4_fit)

  return(mm)
}

gc_fit_moment = function(z,theta_tv_est){ #standardized moments
  fit_tv = likelihood_gc_tv(z,theta_tv_est,type = "n-leverage")
  fit_tv_theta = fit_tv$Theta

  mm_tol = mm_gc_z(fit_tv_theta)

  return(mm_tol)

}

gc_predict = function(fit_variance,data,theta_tv,type = c("n-leverage","leverage","linear")){

  z = fit_variance@fit$z

  t = length(z)

  if(type == "linear"){  #time-varying likelihood
    theta10_tv = theta_tv[1]
    theta11_tv = theta_tv[2]
    theta121_tv = theta_tv[3]
    theta122_tv = theta_tv[3]
    theta13_tv = 0

    theta20_tv = theta_tv[6]
    theta21_tv = theta_tv[7]
    theta221_tv = theta_tv[8]
    theta222_tv = theta_tv[8]
    theta23_tv = 0
  }

  if(type == "leverage"){
    theta10_tv = theta_tv[1]
    theta11_tv = theta_tv[2]
    theta121_tv = theta_tv[3]
    theta122_tv = theta_tv[4]
    theta13_tv = 0

    theta20_tv = theta_tv[6]
    theta21_tv = theta_tv[7]
    theta221_tv = theta_tv[8]
    theta222_tv = theta_tv[9]
    theta23_tv = 0
  }

  if(type == "n-leverage"){
    theta10_tv = theta_tv[1]
    theta11_tv = theta_tv[2]
    theta121_tv = theta_tv[3]
    theta122_tv = theta_tv[4]
    theta13_tv = theta_tv[5]

    theta20_tv = theta_tv[6]
    theta21_tv = theta_tv[7]
    theta221_tv = theta_tv[8]
    theta222_tv = theta_tv[9]
    theta23_tv = theta_tv[10]
  }

  fit_tv = likelihood_gc_tv(z,theta_tv,type = type)
  fit_tv_theta = fit_tv$Theta

  theta1_t1 = tail(fit_tv_theta,1)[1]
  theta2_t1 = tail(fit_tv_theta,1)[2]


  zt = tail(z,1)

  theta1_tp1 = theta10_tv + theta11_tv*theta1_t1 + theta121_tv*(1+ theta13_tv*abs(zt))*max(zt^3,0) + theta122_tv*(1+ theta13_tv*abs(zt))*min(zt^3,0)
  theta2_tp1 = theta20_tv + theta21_tv*theta2_t1 + theta221_tv*(1+ theta23_tv*abs(zt))*max(zt^4,0) + theta222_tv*(1+ theta23_tv*abs(zt))*min(zt^4,0)

  sig_fore = ugarchforecast(fit_variance,data = data,n.ahead = 1)
  sig_tp1 = sig_fore@forecast$sigmaFor
  mu_tp1 = sig_fore@forecast$seriesFor

  skew_ztp1 = theta1_tp1
  kurt_ztp1 = theta2_tp1-3

  mu_xtp1 = mu_tp1
  var_xtp1 = sig_tp1^2
  skew_xtp1 = skew_ztp1*sig_tp1^3
  kurt_xtp1 = kurt_ztp1*sig_tp1^4

  mm_xtp1 = c(mu_xtp1,var_xtp1,skew_xtp1,kurt_xtp1)

  return(mm_xtp1)
}


solve_hess = function(hess){

  id = which(colSums(hess) == 0)
  id_1 = which(colSums(hess) != 0)

  if(length(id) > 0){
    hess_omit = hess[-id,-id]
    std_omit = sqrt(abs(diag(solve(hess_omit))))
    std = vector()
    std[id_1] <- std_omit
    std[id] <- NA
  }

  if(length(id) == 0){
    std = sqrt(abs(diag(solve(hess))))
  }

  return(std)
}


ICA.GARCHSK = function(yy){
  yy = as.matrix(yy)
  n = NCOL(yy)
  t = NROW(yy)
  ica_data = matrix(NA,t,n)
  ica_mean = vector()
  for (bb in 1:n) {
    ica_arma = arima(yy[,bb],order = c(1,0,1))
    ica_data[,bb] <- ica_arma$residuals
    ica_mean[bb] <- predict(ica_arma)$pred
  }
  
  
  jade = ica::icajade(ica_data,nc = n)
  jade_f = jade$S
  jade_B = solve(jade$W)
  
  arcdfit = list()
  arcdfore = list()
  
  f_var_arcd = vector()
  f_skew_arcd = vector()
  f_kurt_arcd = vector()
  
  e_skewpar = e_shapepar = vector()
  for (kk in 1:n){
    e_ik = jade_f[,kk]
    
    spec = rugarch::ugarchspec(mean.model = list(armaOrder = c(0, 0)),variance.model = list(model = 'sGARCH', garchOrder = c(1, 1),variance.targeting = F), distribution = 'sged')
    
    e_fit_variance = rugarch::ugarchfit(spec, e_ik, solver = 'hybrid')
    
    e_z = e_fit_variance@fit$z  #residuals is ht*et, z is et,fitted.value is conditional mean
    
    e_skewpar[kk] = e_fit_variance@fit$coef[4]
    e_shapepar[kk] = e_fit_variance@fit$coef[5]
    
    e_hz = e_fit_variance@fit$residuals
    
    var_e_fit = e_fit_variance@fit$var
    var_e_s = e_hz^2 
    mu_e_fit = e_fit_variance@fit$fitted.values
    
    e_coef_mv = e_fit_variance@fit$coef
    
    e_theta_tv_inl = c(rep(0,5),rep(0,5))
    e_con_tv = optim(par = e_theta_tv_inl,function(theta_tv){con = likelihood_gc_tv_rcpp(e_z,theta_tv,e_theta_est,type = "leverage");return(-con[[1]])}, method = "BFGS", hessian = T)
    
    e_theta_tv_est = e_con_tv$par
    #$residual moment forecasting$
    e_mm_fit = gc_fit_moment(e_z,e_theta_tv_est)
    e_skew_fit = e_mm_fit[,3]
    e_kurt_fit = e_mm_fit[,4]
    
    e_mm_fore = gc_predict(e_fit_variance,e_z,e_theta_tv_est,type = "leverage")
    
    mu_e_fore = e_mm_fore[1]
    var_e_fore = e_mm_fore[2]
    skew_e_fore = e_mm_fore[3]
    kurt_e_fore = e_mm_fore[4]    
    
    f_var_arcd[kk] <-  var_e_fore
    f_skew_arcd[kk] <-  skew_e_fore
    f_kurt_arcd[kk] <-  kurt_e_fore
    
  }
  
  mm_factor_garchsk = list(f_var_arcd,f_skew_arcd,f_kurt_arcd)
  
  return(list(mm_factor_garchsk = mm_factor_garchsk,Beta = jade_B,mu = ica_mean))
  
  
}



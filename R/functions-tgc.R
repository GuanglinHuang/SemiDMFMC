#' @useDynLib SmeiDMFMC
#' @importFrom Rcpp sourceCpp

ugarchspec <- rugarch::ugarchspec
ugarchforecast <- rugarch::ugarchforecast

dtgc = function(x,theta,...){
  theta1 = theta[1]
  theta2 = theta[2]

  gam1 = theta1/sqrt(6)
  gam2 = theta2/sqrt(24)

  lam = 1/(1+gam1^2+gam2^2)

  H3 = (x^3-3*x)/sqrt(6)
  H4 = (x^4-6*x^2+3)/sqrt(24)

  fx = lam*dnorm(x)*(1 + gam1*H3 + gam2*H4)^2

  return(fx)
}

phi = function(x,theta1,theta2,...){
  H3 = (x^3-3*x)/sqrt(6)
  H4 = (x^4-6*x^2+3)/sqrt(24)

  gam1 = theta1/sqrt(6)
  gam2 = theta2/sqrt(24)

  phi = 1+ gam1*H3 + gam2*H4
}

likelihood_tgc_tv = function(z,theta_tv,type = c("n-leverage","leverage","linear"),..){
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

    theta1_t = theta10_tv + theta11_tv*theta1_t1 + theta121_tv*(1+ theta13_tv*abs(zt1))*max(zt1,0) + theta122_tv*(1+ theta13_tv*abs(zt1))*min(zt1,0)
    theta2_t = theta20_tv + theta21_tv*theta2_t1 + theta221_tv*(1+ theta23_tv*abs(zt1))*max(zt1,0) + theta222_tv*(1+ theta23_tv*abs(zt1))*min(zt1,0)

    gam1 = theta1_t/sqrt(6)
    gam2 = theta2_t/sqrt(24)

    lam = 1/(1+gam1^2+gam2^2)

    Ex = 4*lam*gam1*gam2
    Ex2 = 6*lam*gam1^2 + 8*lam*gam2^2 + 1
    Ex3 = 2*sqrt(6)*lam*gam1 + 48*gam1*gam2
    Ex4 = 4*sqrt(6)*lam*gam2 + 72*lam*gam1^2+120*lam*gam2^2+3

    b_theta_t = 1/sqrt(Ex2 - Ex^2)
    a_theta_t = - b_theta_t*Ex

    x_t = (zt - a_theta_t)/b_theta_t

    lt = -log(b_theta_t) + log(lam) - (1/2)*x_t^2 +  2*log(abs(phi(x_t,theta1_t,theta2_t)))

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

likelihood_tgc_tv_rcpp = function(z,theta_tv,theta0,type = c("n-leverage","leverage","linear"),..){
  t = length(z)
  theta10_tv = theta0[1]
  theta20_tv = theta0[2]

  if(type == "linear"){  #time-varying likelihood
    theta_tv[3] <- theta_tv[2]
    theta_tv[7] <- theta_tv[6]
    theta_tv[4] <- 0
    theta_tv[8] <- 0
  }

  if(type == "leverage"){
    theta_tv[4] <- 0
    theta_tv[8] <- 0
  }

  if(type == "n-leverage"){
    theta_tv <- theta_tv
  }

  con = lnf_tgc_tv_tar(z,theta_tv,theta0,length(z))

  return(con)
}

likelihood_tgc_tv_lt_rcpp = function(z,theta_tv,theta0,type = c("n-leverage","leverage","linear"),..){
  t = length(z)
  theta10_tv = theta0[1]
  theta20_tv = theta0[2]

  if(type == "linear"){  #time-varying likelihood
    theta_tv[3] <- theta_tv[2]
    theta_tv[7] <- theta_tv[6]
    theta_tv[4] <- 0
    theta_tv[8] <- 0
  }

  if(type == "leverage"){
    theta_tv[4] <- 0
    theta_tv[8] <- 0
  }

  if(type == "n-leverage"){
    theta_tv <- theta_tv
  }

  con = lnf_tgc_tv_tar_lt(z,theta_tv,theta0,length(z))

  return(con)
}

likelihood_tgc_tv_targeting = function(z,theta_tv,theta0,type = c("n-leverage","leverage","linear"),..){
  t = length(z)
  theta10_tv = theta0[1]
  theta20_tv = theta0[2]

  if(type == "linear"){  #time-varying likelihood

    theta11_tv = theta_tv[1]
    theta121_tv = theta_tv[2]
    theta122_tv = theta121_tv
    theta13_tv = 0


    theta21_tv = theta_tv[5]
    theta221_tv = theta_tv[6]
    theta222_tv = theta221_tv
    theta23_tv = 0
  }

  if(type == "leverage"){
    theta11_tv = theta_tv[1]
    theta121_tv = theta_tv[2]
    theta122_tv = theta_tv[3]
    theta13_tv = 0

    theta21_tv = theta_tv[5]
    theta221_tv = theta_tv[6]
    theta222_tv = theta_tv[7]
    theta23_tv = 0
  }

  if(type == "n-leverage"){
    theta11_tv = theta_tv[1]
    theta121_tv = theta_tv[2]
    theta122_tv = theta_tv[3]
    theta13_tv = theta_tv[4]

    theta21_tv = theta_tv[5]
    theta221_tv = theta_tv[6]
    theta222_tv = theta_tv[7]
    theta23_tv = theta_tv[8]
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

  for (tt in 2:t){
    zt1 = z[tt-1]
    zt = z[tt]

    # if(zt1 > 0){
    #   theta1_t = theta10_tv + theta11_tv*theta1_t1 + theta121_tv*(1 + theta13_tv*abs(zt1))*zt1;
    #   theta2_t = theta20_tv + theta21_tv*theta2_t1 + theta221_tv*(1 + theta23_tv*abs(zt1))*zt1;
    # }
    # if(zt1 <= 0){
    #   theta1_t = theta10_tv + theta11_tv*theta1_t1 + theta122_tv*(1 + theta13_tv*abs(zt1))*zt1;
    #   theta2_t = theta20_tv + theta21_tv*theta2_t1 + theta222_tv*(1 + theta23_tv*abs(zt1))*zt1;
    # }

    theta1_t = theta10_tv + theta11_tv*theta1_t1 + theta121_tv*(1+ theta13_tv*abs(zt1))*max(zt1,0) + theta122_tv*(1+ theta13_tv*abs(zt1))*min(zt1,0)
    theta2_t = theta20_tv + theta21_tv*theta2_t1 + theta221_tv*(1+ theta23_tv*abs(zt1))*max(zt1,0) + theta222_tv*(1+ theta23_tv*abs(zt1))*min(zt1,0)

    gam1 = theta1_t/sqrt(6)
    gam2 = theta2_t/sqrt(24)

    lam = 1/(1+gam1^2+gam2^2)

    Ex = 4*lam*gam1*gam2
    Ex2 = 6*lam*gam1^2 + 8*lam*gam2^2 + 1
    Ex3 = 2*sqrt(6)*lam*gam1 + 48*gam1*gam2
    Ex4 = 4*sqrt(6)*lam*gam2 + 72*lam*gam1^2+120*lam*gam2^2+3

    b_theta_t = 1/sqrt(Ex2 - Ex^2)
    a_theta_t = - b_theta_t*Ex

    x_t = (zt - a_theta_t)/b_theta_t

    H3 = (x_t*x_t*x_t-3*x_t)/sqrt(6)
    H4 = (x_t*x_t*x_t*x_t-6*x_t*x_t+3)/sqrt(24)
    phi = 1 + gam1*H3 + gam2*H4

    lt = -log(b_theta_t) + log(lam) - (1/2)*x_t^2 +  2*log(abs(phi))

    lt_tol = c(lt_tol,lt)

    theta1_tv[tt] = theta1_t
    theta2_tv[tt] = theta2_t

    theta1_t1 = theta1_t
    theta2_t1 = theta2_t
  }

  LL = sum(lt_tol)

  con = list(Lnf = LL,Theta = cbind(theta1_tv,theta2_tv),lt_tol)

  return(con)
}

likelihood_tgc_tv_targeting_rcpp = function(z,theta_tv,theta0,type = c("n-leverage","leverage","linear"),..){
  t = length(z)
  theta10_tv = theta0[1]
  theta20_tv = theta0[2]

  if(type == "linear"){  #time-varying likelihood
    theta_tv[3] <- theta_tv[2]
    theta_tv[7] <- theta_tv[6]
    theta_tv[4] <- 0
    theta_tv[8] <- 0
  }

  if(type == "leverage"){
    theta_tv[4] <- 0
    theta_tv[8] <- 0
  }

  if(type == "n-leverage"){
    theta_tv <- theta_tv
  }

  con = lnf_tgc_tv_tar(z,theta_tv,theta0,length(z))

  return(con)
}

likelihood_tgc_tv_targeting_lt_rcpp = function(z,theta_tv,theta0,type = c("n-leverage","leverage","linear"),..){
  t = length(z)
  theta10_tv = theta0[1]
  theta20_tv = theta0[2]

  if(type == "linear"){  #time-varying likelihood
    theta_tv[3] <- theta_tv[2]
    theta_tv[7] <- theta_tv[6]
    theta_tv[4] <- 0
    theta_tv[8] <- 0
  }

  if(type == "leverage"){
    theta_tv[4] <- 0
    theta_tv[8] <- 0
  }

  if(type == "n-leverage"){
    theta_tv <- theta_tv
  }

  con = lnf_tgc_tv_tar_lt(z,theta_tv,theta0,length(z))

  return(con)
}

likelihood_tgc_con = function(z,theta,..){  #constant likehihood

  t = length(z)

  theta1 = theta[1]
  theta2 = theta[2]

  gam1 = theta1/sqrt(6)
  gam2 = theta2/sqrt(24)

  lam = 1/(1+gam1^2+gam2^2)

  Ex = 4*lam*gam1*gam2
  Ex2 = 6*lam*gam1^2 + 8*lam*gam2^2 + 1
  Ex3 = 2*sqrt(6)*lam*gam1 + 48*gam1*gam2
  Ex4 = 4*sqrt(6)*lam*gam2 + 72*lam*gam1^2+120*lam*gam2^2+3

  b_theta = 1/sqrt(Ex2 - Ex^2)
  a_theta = - b_theta*Ex

  x = (z - a_theta)/b_theta

  lt = -log(b_theta) + log(lam) - (1/2)*x^2 +  2*log(abs(phi(x,theta1,theta2)))

  LL = sum(lt)

  con = list(Lnf = LL,Theta = c(theta1,theta2),lt)

  return(con)
}

#tgc moments

mm_tgc_z = function(theta){ #standardized moments
  theta1 = theta[,1]
  theta2 = theta[,2]
  gam1 = theta1/sqrt(6)
  gam2 = theta2/sqrt(24)
  lam = 1/(1 + gam1^2 + gam2^2)

  mm1 = 4*lam*gam1*gam2
  mm2 = 6*lam*gam1^2 + 8*lam*gam2^2 + 1
  mm3 = 2*sqrt(6)*lam*gam1 + 48*lam*gam1*gam2
  mm4 = 4*sqrt(6)*lam*gam2 + 72*lam*gam1^2 + 120*lam*gam2^2 + 3

  b = 1/sqrt(mm2 - mm1^2)
  a = - b*mm1

  m1_fit = 0
  m2_fit = 1
  m3_fit = a^3+3*a^2*b*mm1 + 3*a*b^2*mm2+b^3*mm3
  m4_fit = a^4+4*a^3*b*mm1 + 6*a^2*b^2*mm2 + 4*a*b^3*mm3+b^4*mm4 - 3

  mm = cbind(m1_fit,m2_fit,m3_fit,m4_fit)

  return(mm)
}



tgc_fit_moment = function(z,theta_tv_est){
  fit_tv = likelihood_tgc_tv(z,theta_tv_est,type = "n-leverage")
  fit_tv_theta = fit_tv$Theta

  mm_tol = mm_tgc_z(fit_tv_theta)

  return(mm_tol)

}

tgc_predict = function(fit_variance,data,theta_tv,type = c("n-leverage","leverage","linear")){

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

  fit_tv = likelihood_tgc_tv(z,theta_tv,type = type)
  fit_tv_theta = fit_tv$Theta

  theta1_t1 = tail(fit_tv_theta,1)[1]
  theta2_t1 = tail(fit_tv_theta,1)[2]


  zt = tail(z,1)

  theta1_tp1 = theta10_tv + theta11_tv*theta1_t1 + theta121_tv*(1+ theta13_tv*abs(zt))*max(zt,0) + theta122_tv*(1+ theta13_tv*abs(zt))*min(zt,0)
  theta2_tp1 = theta20_tv + theta21_tv*theta2_t1 + theta221_tv*(1+ theta23_tv*abs(zt))*max(zt,0) + theta222_tv*(1+ theta23_tv*abs(zt))*min(zt,0)

  gam1 = theta1_tp1/sqrt(6)
  gam2 = theta2_tp1/sqrt(24)

  lam = 1/(1+gam1^2+gam2^2)


  Ez1 = 4*lam*gam1*gam2
  Ez2 = 6*lam*gam1^2 + 8*lam*gam2^2 + 1
  Ez3 = 2*sqrt(6)*lam*gam1 + 48*gam1*gam2
  Ez4 = 4*sqrt(6)*lam*gam2 + 72*lam*gam1^2+120*lam*gam2^2 + 3

  #skew_ztp1 = Ez3
  #kurt_ztp1 = Ez4 - 3

  sig_fore = ugarchforecast(fit_variance,data = data,n.ahead = 1)
   sig_tp1 = sig_fore@forecast$sigmaFor
    mu_tp1 = sig_fore@forecast$seriesFor

  skew_ztp1 = (Ez3 - 3*Ez1*Ez2 + 2*(Ez1)^3)/(Ez2-Ez1^2)^1.5
  kurt_ztp1 = ((Ez4 - 4*Ez1*Ez3 + 6*Ez2*Ez1^2 - 3*Ez1^4)/(Ez2-Ez1^2)^2-3)

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
    std_omit = sqrt(abs(diag(MASS::ginv(hess_omit))))
    std = vector()
    std[id_1] <- std_omit
    std[id] <- NA
  }

  if(length(id) == 0){
    std = sqrt(abs(diag(MASS::ginv(hess))))
  }

  return(std)
}

TVSNP.test = function(z,result_tv,result_con,type = "leverage"){
  
  theta_est = result_con$snp.cof
  lt_csnp =  likelihood_tgc_con(z,theta_est)[[3]][-1]
  
  theta_tv_est = result_tv$tgc.cof
  type = type
  lt_tvsnp = c(likelihood_tgc_tv_lt_rcpp(z,theta_tv_est[c(2:5,7:10)],theta0 = theta_tv_est[c(1,6)],type = type))
  
  LLC = sum(lt_csnp) 
  LLTV = sum(lt_tvsnp)
  
  wt = sqrt(mean((lt_tvsnp-lt_csnp)^2)-(mean(lt_tvsnp-lt_csnp))^2)
  LR_f = (LLTV-LLC)/(sqrt(length(lt_csnp))*wt)
  
  
  return(list(LR_f = LR_f,LLTV = LLTV,LLC = LLC))
}

TGC_est = function(ff,var.model = 'sGARCH',var.targeting = F,var.distribution = 'sged',
                   tgc.type  = "leverage",tgc.targeting = F, mean.model = list(armaOrder = c(1,1)),CTGC = FALSE,rep_sim = 10){
  
  type  = tgc.type
  
  spec = rugarch::ugarchspec(mean.model = mean.model,variance.model = list(model = var.model, garchOrder = c(1, 1),variance.targeting = var.targeting), distribution = var.distribution)
  
  fit_variance = rugarch::ugarchfit(spec, ff, solver = 'hybrid')
  
  z = fit_variance@fit$z  #residuals is ht*et, z is et,fitted.value is conditional mean
  
  hz = fit_variance@fit$residuals
  
  var_x_fit = fit_variance@fit$var
  var_x_s = hz^2 
  mu_x_fit = fit_variance@fit$fitted.values
  
  coef_mv = fit_variance@fit$coef
  
  if(CTGC == T){
    #$constant theta estimation$

      theta_con_inl = runif(2,-1,1)

      constant_est = Rsolnp::gosolnp(theta_con_inl,fun = function(theta){
        con = likelihood_tgc_con(z,theta)
        return(-con[[1]])
      },LB = rep(-2,2),UB = rep(2,2), n.sim = 5000,n.restarts = rep_sim, control = list(trace = 0))
  
      
    std = sqrt(diag(MASS::ginv(constant_est$hessian)))
    t_stat = constant_est$pars/std
    p_value = 2*(1 - pnorm(abs(t_stat)))
    
    constant_con = cbind(constant_est$pars,std,t_stat,p_value) #constant parameters estimation #return
    row.names(constant_con) <- c("theta1","theta2")
    theta_est = constant_con[,1] #return
    
    #moment forecasting
    theta_tv_inl = c(theta_est[1],rep(0,4),theta_est[2],rep(0,4))
    
    result_moment = tgc_predict(fit_variance,ff,theta_tv_inl,type = type)
    
    Lnf_con = likelihood_tgc_con(z,theta_est)$Lnf
    
    result_con = list(
      snp.cof = theta_est, constant_con = constant_con, result_moment = list(mm.fore = result_moment),
      Lnf_con = Lnf_con
    )
    
    return(result_con)
  }
  
  if(CTGC == F){
    
    theta_con_inl = runif(2,-1,1)
    
    constant_est = Rsolnp::solnp(theta_con_inl,fun = function(theta){
      con = likelihood_tgc_con(z,theta)
      return(-con[[1]])
    },LB = rep(-1,2),UB = rep(1,2), control = list(trace = 0))
    
    
    std = sqrt(diag(MASS::ginv(constant_est$hessian)))
    t_stat = constant_est$pars/std
    p_value = 2*(1 - pnorm(abs(t_stat)))
    
    constant_con = cbind(constant_est$pars,std,t_stat,p_value) #constant parameters estimation #return
    row.names(constant_con) <- c("theta1","theta2")
    theta_est = constant_con[,1] #return
    
    #moment forecasting
    theta_tv_inl = c(theta_est[1],rep(0,4),theta_est[2],rep(0,4))
    
    theta_tv_inl_con = theta_tv_inl
    
    result_moment = tgc_predict(fit_variance,ff,theta_tv_inl,type = type)
    
    Lnf_con = likelihood_tgc_con(z,theta_est)$Lnf
    
    result_con = list(
      snp.cof = theta_est, constant_con = constant_con, result_moment = list(mm.fore = result_moment),
      Lnf_con = Lnf_con
    )
    
    #$time varying theta estimation$
    if(!tgc.targeting){
      
      if(tgc.type == "linear"){
        ar1_inl = runif(1,0.5,0.9)
        ar2_inl = runif(1,0.5,0.9)
        theta_tv_inl = c(0,ar1_inl,0,0,0,ar2_inl,0,0)
        
        LB = rep(-1,8)
        UB = rep(1,8)
        
        con_tv = Rsolnp::gosolnp(pars = theta_tv_inl,fixed = c(4,8),function(sss){
          theta_tv = sss[c(2,3,3,4,6,7,7,8)];
          theta0   = c(sss[1],sss[5]);
          con      = likelihood_tgc_tv_rcpp(z,theta_tv,theta0,type = "linear"); 
          return(-con)},LB = LB, UB = UB,n.sim = 20000,n.restarts = rep_sim, control = list(trace = 0))
        
        std_tv = solve_hess(con_tv$hessian)
        t_stat_tv = con_tv$pars/std_tv
        p_value_tv = 2*(1-pnorm(abs(t_stat_tv)))
        tv_con = cbind(con_tv$pars,std_tv,t_stat_tv,p_value_tv)#tv parameters estimation
        row.names(tv_con)<-c("omega1","alpha1","beta1","gamma1","omega2","alpha2","beta2","gamma2")
        
        theta_tv_est = con_tv$pars[c(1,2,3,3,4,5,6,6,7,8)]
      }
      
      if(tgc.type == "leverage"){
        
        ar1_inl = runif(1,0.5,0.9)
        ar2_inl = runif(1,0.5,0.9)
        theta_tv_inl = c(0,ar1_inl,0,0,0,0,ar2_inl,0,0,0)
        
        LB = c(-0.5,-1,-1,-1,-0.3,-0.5,-1,-1,-1,-0.3)
        UB = c(0.5,1,1,1,0.3,0.5,1,1,1,0.3)
        
        con_tv = Rsolnp::gosolnp(pars = theta_tv_inl,fixed = c(5,10),function(sss){
          theta_tv = sss[-c(1,6)];
          theta0   = c(sss[1],sss[6]);
          con      = likelihood_tgc_tv_rcpp(z,theta_tv,theta0,type = "leverage"); 
          return(-con)},LB = LB, UB = UB,n.sim = 20000,n.restarts = rep_sim, control = list(trace = 0))
        
        std_tv = solve_hess(con_tv$hessian)
        t_stat_tv = con_tv$pars/std_tv
        p_value_tv = 2*(1-pnorm(abs(t_stat_tv)))
        tv_con = cbind(con_tv$pars,std_tv,t_stat_tv,p_value_tv)#tv parameters estimation
        row.names(tv_con)<-c("omega1","alpha1","beta11","beta12","gamma1","omega2","alpha2","beta21","beta22","gamma2")
        
        theta_tv_est = con_tv$pars
      }
      
      if(tgc.type == "n-leverage"){
        ar1_inl = runif(1,0.5,0.8)
        ar2_inl = runif(1,0.5,0.8)
        theta_tv_inl = c(0,ar1_inl,0,0,0,0,ar2_inl,0,0,0)
        
        LB = c(-0.5,-1,-1,-1,-0.3,-0.5,-1,-1,-1,-0.3)
        UB = c(0.5,1,1,1,0.3,0.5,1,1,1,0.3)
        
        con_tv = Rsolnp::gosolnp(pars = theta_tv_inl,function(sss){
          theta_tv = sss[-c(1,6)];
          theta0   = c(sss[1],sss[6]);
          con      = likelihood_tgc_tv_rcpp(z,theta_tv,theta0,type = "n-leverage"); 
          return(-con)},LB = LB, UB = UB,n.sim = 20000,n.restarts = rep_sim,control = list(trace = 0))
        
        std_tv = solve_hess(con_tv$hessian)
        t_stat_tv = con_tv$pars/std_tv
        p_value_tv = 2*(1-pnorm(abs(t_stat_tv)))
        tv_con = cbind(con_tv$pars,std_tv,t_stat_tv,p_value_tv)#tv parameters estimation
        row.names(tv_con)<-c("omega1","alpha1","beta11","beta12","gamma1","omega2","alpha2","beta21","beta22","gamma2")
        
        theta_tv_est = con_tv$pars
      }
 
    }
    
    # targeting estimation
    if(tgc.targeting){
      if(tgc.type == "linear"){
        ar1_inl = runif(1,0.5,0.8)
        ar2_inl = runif(1,0.5,0.8)
        
        theta_tv_inl = c(theta_tv_inl_con[1],ar1_inl,0,0,theta_tv_inl_con[6],ar2_inl,0,0)
        
        LB = c(-0.5,-1,-1,-1,-0.3,-0.5,-1,-1,-1,-0.3)
        UB = c(0.5,1,1,1,0.3,0.5,1,1,1,0.3)
        
        con_tv = Rsolnp::gosolnp(pars = theta_tv_inl,fixed = c(1,4,5,8),function(sss){
          theta_tv = sss[c(2,3,3,4,6,7,7,8)];
          theta0   = c(theta_tv_inl_con[1],theta_tv_inl_con[6]);
          con      = likelihood_tgc_tv_targeting_rcpp(z,theta_tv,theta0,type = "linear"); 
          return(-con)},LB = LB, UB = UB,n.sim = 20000,n.restarts = rep_sim,control = list(trace = 0))
        
        std_tv = solve_hess(con_tv$hessian)
        t_stat_tv = con_tv$pars/std_tv
        p_value_tv = 2*(1-pnorm(abs(t_stat_tv)))
        tv_con = cbind(con_tv$pars,std_tv,t_stat_tv,p_value_tv)#tv parameters estimation
        tv_con[1,] = constant_con[1,]
        tv_con[5,] = constant_con[2,]
        row.names(tv_con) <- c("theta1","alpha1","beta1","gamma1","theta2","alpha2","beta2","gamma2")
        
        theta_tv_est = con_tv$pars[c(1,2,3,3,4,5,6,6,7,8)]
        names(theta_tv_est) = c("theta1","alpha1","beta11","beta12","gamma1","theta2","alpha2","beta21","beta22","gamma2")
        
      }
      
      if(tgc.type == "leverage"){
        ar1_inl = runif(1,0.5,0.8)
        ar2_inl = runif(1,0.5,0.8)
        theta_tv_inl = c(theta_tv_inl_con[1],ar1_inl,0,0,0,theta_tv_inl_con[6],ar2_inl,0,0,0)
        
        LB = c(-0.5,-1,-1,-1,-0.3,-0.5,-1,-1,-1,-0.3)
        UB = c(0.5,1,1,1,0.3,0.5,1,1,1,0.3)
        
        con_tv = Rsolnp::gosolnp(pars = theta_tv_inl,fixed = c(1,5,6,10),function(sss){
          theta_tv = sss[-c(1,6)];
          theta0   = c(theta_tv_inl_con[1],theta_tv_inl_con[6]);
          con      = likelihood_tgc_tv_rcpp(z,theta_tv,theta0,type = "leverage"); 
          return(-con)},LB = LB, UB = UB,n.sim = 20000,n.restarts = rep_sim,control = list(trace = 0))
        
        std_tv = solve_hess(con_tv$hessian)
        t_stat_tv = con_tv$pars/std_tv
        p_value_tv = 2*(1-pnorm(abs(t_stat_tv)))
        tv_con = cbind(con_tv$pars,std_tv,t_stat_tv,p_value_tv)#tv parameters estimation
        tv_con[1,] = constant_con[1,]
        tv_con[6,] = constant_con[2,]
        row.names(tv_con)<-c("theta1","alpha1","beta11","beta12","gamma1","theta2","alpha2","beta21","beta22","gamma2")
        
        theta_tv_est = con_tv$pars
        names(theta_tv_est) = c("theta1","alpha1","beta11","beta12","gamma1","theta2","alpha2","beta21","beta22","gamma2")
        
      }
      
      if(tgc.type == "n-leverage"){
        ar1_inl = runif(1,0.5,0.8)
        ar2_inl = runif(1,0.5,0.8)
        theta_tv_inl = c(0,ar1_inl,0,0,0,0,ar2_inl,0,0,0)
        
        LB = c(-0.5,-1,-1,-1,-0.3,-0.5,-1,-1,-1,-0.3)
        UB = c(0.5,1,1,1,0.3,0.5,1,1,1,0.3)
        
        con_tv = Rsolnp::gosolnp(pars = theta_tv_inl,fixed = c(1,6),function(sss){
          theta_tv = sss[-c(1,6)];
          theta0   = c(theta_tv_inl_con[1],theta_tv_inl_con[6]);
          con      = likelihood_tgc_tv_rcpp(z,theta_tv,theta0,type = "n-leverage"); 
          return(-con)},LB = LB, UB = UB,n.sim = 20000,n.restarts = rep_sim,control = list(trace = 0))
        
        std_tv = solve_hess(con_tv$hessian)
        t_stat_tv = con_tv$pars/std_tv
        p_value_tv = 2*(1-pnorm(abs(t_stat_tv)))
        tv_con = cbind(con_tv$pars,std_tv,t_stat_tv,p_value_tv)#tv parameters estimation
        tv_con[1,] = constant_con[1,]
        tv_con[6,] = constant_con[2,]
        row.names(tv_con)<-c("theta1","alpha1","beta11","beta12","gamma1","theta2","alpha2","beta21","beta22","gamma2")
        
        theta_tv_est = con_tv$pars
        names(theta_tv_est) = c("theta1","alpha1","beta11","beta12","gamma1","theta2","alpha2","beta21","beta22","gamma2")
        
      }    
}
    
    
    result_tv = list(
      tgc.cof = theta_tv_est, tv_con = tv_con
      )
    
    #$moment forecasting$
    mm_fit = tgc_fit_moment(z,theta_tv_est) #return
    
    skew_fit = mm_fit[,3]
    kurt_fit = mm_fit[,4]
    
    mm_fore = tgc_predict(fit_variance,ff,theta_tv_est,type = type)
    
    #likelihood 
    Lnf_tv = likelihood_tgc_tv_rcpp(z,theta_tv_est[-c(1,6)],theta_tv_est[c(1,6)],type = type)
    
    result_moment = list(
      mm.fit = mm_fit,
      mm.fore = mm_fore
    )
    
    return(list(result_con = result_con,
                result_tv = result_tv,
                result_moment = result_moment,
                fit_variance =fit_variance,
                Lnf_tv = Lnf_tv)
    )
  }
  
}

MFTGC_est = function(ff,var.model = 'sGARCH',var.targeting = F,var.distribution = 'sged', VAR = T,Corr.Struture = c("ica","dcc","copula"), dcc.model = "DCC", 
                     copula.model = list(copula = "mvt", method = "ML", time.varying = T, transformation = "spd"),
                     tgc.type  = "leverage",tgc.targeting = F,mean.model = list(armaOrder = c(0, 0)),CTGC = FALSE,rep_sim = 10){
  
  result_factors_ica = NULL
  result_factors_dcc = NULL
  result_factors_copula = NULL
  
  ica_moments = NULL
  dcc_moments = NULL
  copula_moments = NULL
  
  q = NCOL(ff)
  ff = as.matrix(ff)
  
  if("ica" %in% Corr.Struture){
    
    if(VAR == T){
      gogarch = rmgarch::gogarchspec(mean.model = list(model = "VAR", robust = FALSE, 
                                                       lag = 1, lag.max = NULL, lag.criterion = "FPE", 
                                                       external.regressors = NULL, 
                                                       robust.control = list("gamma" = 0.25, "delta" = 0.01, "nc" = 10, "ns" = 500)), 
                                     variance.model = list(model = var.model, garchOrder = c(1,1), submodel = NULL, 
                                                           variance.targeting = var.targeting), distribution.model = "manig", 
                                     ica = "fastica", ica.fix = list(A = NULL, K = NULL)) 
    }else{
      gogarch = rmgarch::gogarchspec(mean.model = list(model = "constant", robust = FALSE, 
                                                       lag = 1, lag.max = NULL, lag.criterion = "FPE", 
                                                       external.regressors = NULL, 
                                                       robust.control = list("gamma" = 0.25, "delta" = 0.01, "nc" = 10, "ns" = 500)), 
                                     variance.model = list(model = var.model, garchOrder = c(1,1), submodel = NULL, 
                                                           variance.targeting = var.targeting), distribution.model = "manig", 
                                     ica = "fastica", ica.fix = list(A = NULL, K = NULL))  
    }

    fit_gogarch = rmgarch::gogarchfit(gogarch,ff,out.sample = 0, solver = "solnp", gfun = "tanh")
    
    fore_gogarch = rmgarch::gogarchforecast(fit_gogarch,n.ahead = 1)

    ica_B = fit_gogarch@mfit$A

    lf = fit_gogarch@mfit$residuals%*%t(solve(ica_B))
    
    #independent component moment estimation
    var_mf_fore = matrix(0,q,q)
    skew_mf_fore = matrix(0,q,q^2)
    kurt_mf_fore = matrix(0,q,q^3)
    
    lnf_lf = vector()
    lnf_lf_con = vector()
    
    con_factor <- list()
    
    for (gg in 1:q) {
      lf_snp = TGC_est(lf[,gg],mean.model = list(armaOrder = c(0, 0)),tgc.type = tgc.type,CTGC = CTGC,tgc.targeting = tgc.targeting,rep_sim = rep_sim)
      lf_snp$Lnf_tv
      lf_snp$result_con$Lnf_con
      
      con_factor[[gg]] <- lf_snp
      
      mu_factor_fore = lf_snp$result_moment$mm.fore[1]
      var_factor_fore = lf_snp$result_moment$mm.fore[2]
      skew_factor_fore = lf_snp$result_moment$mm.fore[3]
      kurt_factor_fore = lf_snp$result_moment$mm.fore[4]
      
      lnf_lf[gg] = lf_snp$Lnf_tv
      lnf_lf_con[gg] = lf_snp$result_con$Lnf_con
      
      var_mf_fore[gg,gg] = var_factor_fore
      skew_mf_fore[gg,gg+(gg-1)*q] = skew_factor_fore
      kurt_mf_fore[gg,gg+(gg-1)*q+(gg-1)*q^2] = kurt_factor_fore
    }
    
    var_f_fore = ica_B%*%var_mf_fore%*%t(ica_B)
    skew_f_fore = ica_B%*%skew_mf_fore%*%(t(ica_B)%x%t(ica_B))
    kurt_f_fore = ica_B%*%kurt_mf_fore%*%(t(ica_B)%x%t(ica_B)%x%t(ica_B))
    
    result_factors_ica = list(result_mgarch = fit_gogarch,fore_mgarch = fore_gogarch,result_snp = con_factor,factor_moments = list(var_f_fore,skew_f_fore,kurt_f_fore))
    
    ica_moments = list(var_f_fore,skew_f_fore,kurt_f_fore)
  }
  
  if("dcc" %in% Corr.Struture){
  
    uspec = rugarch::ugarchspec(mean.model = mean.model,
                                variance.model = list(model = var.model, garchOrder = c(1, 1),
                                                      variance.targeting = var.targeting), distribution = var.distribution) #univariate spec
    
    mspec = rugarch::multispec( replicate(q, uspec) )
    
    if(VAR == T){
      dccgarch = rmgarch::dccspec(uspec = mspec, VAR = T, model = dcc.model,distribution = "mvt")
    }else{
      dccgarch = rmgarch::dccspec(uspec = mspec, VAR = F, model = dcc.model,distribution = "mvt")
    }
    
    
    fit_dccgarch = rmgarch::dccfit(dccgarch,ff)
    
    fore_dccgarch = rmgarch::dccforecast(fit_dccgarch,n.ahead = 1)
    
    var_f_fore = as.matrix(fore_dccgarch@mforecast$H[[1]][,,1])
    
    lf = fit_dccgarch@mfit$stdresid
    
    #independent component moment estimation
    q = NCOL(lf)
    var_mf_fore = matrix(0,q,q)
    skew_mf_fore = matrix(0,q,q^2)
    kurt_mf_fore = matrix(0,q,q^3)
    
    lnf_lf = vector()
    lnf_lf_con = vector()
    
    con_factor <- list()
    
    for (gg in 1:NCOL(ff)) {
      lf_snp = TGC_est(lf[,gg],mean.model = list(armaOrder = c(0, 0)),tgc.type = tgc.type,tgc.targeting = tgc.targeting,CTGC = CTGC,rep_sim = rep_sim)
      
      con_factor[[gg]] <- lf_snp
      
      mu_factor_fore = lf_snp$result_moment$mm.fore[1]
      var_factor_fore = lf_snp$result_moment$mm.fore[2]
      skew_factor_fore = lf_snp$result_moment$mm.fore[3]
      kurt_factor_fore = lf_snp$result_moment$mm.fore[4]
      
      lnf_lf[gg] = lf_snp$Lnf_tv
      lnf_lf_con[gg] = lf_snp$result_con$Lnf_con
      
      var_mf_fore[gg,gg] = var_factor_fore/var_factor_fore
      skew_mf_fore[gg,gg+(gg-1)*q] = skew_factor_fore/(var_factor_fore)^1.5
      kurt_mf_fore[gg,gg+(gg-1)*q+(gg-1)*q^2] = kurt_factor_fore/(var_factor_fore)^2
    }

    sig_f_fore = Mat.k(var_f_fore,1/2)
    
    var_f_fore = sig_f_fore%*%var_mf_fore%*%sig_f_fore
    skew_f_fore = sig_f_fore%*%skew_mf_fore%*%(sig_f_fore%x%sig_f_fore)
    kurt_f_fore = sig_f_fore%*%kurt_mf_fore%*%(sig_f_fore%x%sig_f_fore%x%sig_f_fore)
    
    result_factors_dcc = list(result_mgarch = fit_dccgarch,fore_mgarch = fore_dccgarch,result_snp = con_factor,factor_moments = list(var_f_fore,skew_f_fore,kurt_f_fore))
    dcc_moments = list(var_f_fore,skew_f_fore,kurt_f_fore)
  }
  
  if("copula" %in% Corr.Struture){
    
    uspec = rugarch::ugarchspec(mean.model = mean.model,
                                variance.model = list(model = var.model, garchOrder = c(1, 1),
                                                      variance.targeting = var.targeting), distribution = var.distribution) #univariate spec
    
    mspec = rugarch::multispec( replicate(q, uspec) )
    
    if(VAR == T){
      copulagarch = rmgarch::cgarchspec(uspec = mspec, VAR = T, distribution.model = copula.model)
    }else{
      copulagarch = rmgarch::cgarchspec(uspec = mspec, VAR = F, distribution.model = copula.model)
    }
    
    
    fit_copulagarch = rmgarch::cgarchfit(copulagarch,ff)
    
    # use simulation to forecast the Ht
    sim2 = rmgarch::cgarchsim(fit_copulagarch, n.sim = 1, m.sim = 1000, startMethod = "sample")
    
    mu_f_fore = 0
    for (ij in 1:1000) {
      mu_f_fore = mu_f_fore + sim2@msim$simX[[ij]]
    }
    mu_f_fore = mu_f_fore/1000
    
    var_f_fore = sim2@msim$simH[[1]][,,1]
    
    lf = fit_copulagarch@mfit$stdresid
    
    #independent component moment estimation
  
    var_mf_fore = matrix(0,q,q)
    skew_mf_fore = matrix(0,q,q^2)
    kurt_mf_fore = matrix(0,q,q^3)
    
    lnf_lf = vector()
    lnf_lf_con = vector()
    
    con_factor <- list()
    
    for (gg in 1:NCOL(ff)) {
      lf_snp = TGC_est(lf[,gg],mean.model = list(armaOrder = c(0, 0)),tgc.type = tgc.type,CTGC = CTGC,tgc.targeting = tgc.targeting,rep_sim = rep_sim)
      
      con_factor[[gg]] <- lf_snp
      
      mu_factor_fore = lf_snp$result_moment$mm.fore[1]
      var_factor_fore = lf_snp$result_moment$mm.fore[2]
      skew_factor_fore = lf_snp$result_moment$mm.fore[3]
      kurt_factor_fore = lf_snp$result_moment$mm.fore[4]
      
      lnf_lf[gg] = lf_snp$Lnf_tv
      lnf_lf_con[gg] = lf_snp$result_con$Lnf_con
      
      var_mf_fore[gg,gg] = var_factor_fore/var_factor_fore
      skew_mf_fore[gg,gg+(gg-1)*q] = skew_factor_fore/(var_factor_fore)^1.5
      kurt_mf_fore[gg,gg+(gg-1)*q+(gg-1)*q^2] = kurt_factor_fore/(var_factor_fore)^2
    }
    
    
    sig_f_fore = Mat.k(var_f_fore,1/2)
    
    var_f_fore = sig_f_fore%*%var_mf_fore%*%sig_f_fore
    skew_f_fore = sig_f_fore%*%skew_mf_fore%*%(sig_f_fore%x%sig_f_fore)
    kurt_f_fore = sig_f_fore%*%kurt_mf_fore%*%(sig_f_fore%x%sig_f_fore%x%sig_f_fore)
  
    result_factors_copula = list(result_mgarch = fit_copulagarch,fore_mgarch = list(mu = mu_f_fore), result_snp = con_factor,factor_moments = list(var_f_fore,skew_f_fore,kurt_f_fore))
    
    copula_moments = list(var_f_fore,skew_f_fore,kurt_f_fore)
  }
  
  
  return(list(result_factors_ica      =   result_factors_ica,
              result_factors_dcc      =   result_factors_dcc,
              result_factors_copula   =   result_factors_copula,
              factor_moments          =   list(ica_moments    = ica_moments,
                                               dcc_moments    = dcc_moments,
                                               copula_moments = copula_moments)))
  
}
  
  



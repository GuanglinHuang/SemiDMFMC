#' @useDynLib SmeiDMFMC
#' @importFrom Rcpp sourceCpp


npscoefbw <- np::npscoefbw

npest_tv = function(u,h,y,x,z){
  XX = vector()
  W = vector()
  t = length(z)
  YY = y[2:t,]
  FF = x
  ZZ = z
  XX = cbind(rep(1,t-1),FF[2:t,],ZZ[1:(t-1)]-u,FF[2:t,]*(ZZ[1:(t-1)]-u))
  WW = apply(matrix(ZZ[1:(t-1)],t-1,1),1,function(x){K((u-x)/h)/h})
  WW = diag(WW,t-1,t-1)
  EST = MASS::ginv(t(XX)%*%WW%*%XX)%*%t(XX)%*%WW%*%YY
  return(EST)
}


cv_h1_tv <- function(M = 0.1,X,FF,kmax = 100,trimp = 0.2){
  cv = vector()
  bd_t_t = vector()
  t = length(FF)
  MM = ceiling(M*t)-1
  w  =  abs(FF[1:(t-1)] - FF[t])
  if(length(which(duplicated(w) == T))>0){
    w  = w[-which(duplicated(w) == T)]
  }
  w = sort(w)
  wl = quantile(w,trimp)
  wu = quantile(w,1-trimp)
  ww = seq(from = wl,to = wu,length.out  = kmax+1)
  for (k in 1:kmax){
    sq_total = 0
    bd = ww[k]
    for (tt in (t-MM):t){
      a = npest_tv(FF[tt-1],bd,as.matrix(X[1:(tt-1),]),as.matrix(FF[1:(tt-1)]))
      cof = a[,1:2]
      sq = fnorm(X[tt,]- cof%*%c(1,FF[tt]))
      sq_total = sq_total + sq
    }
    bd_t_t[k] <- bd
    cv[k] <- sq_total
    
  }
  if(T %in% is.nan(cv)){
    bd_t_t <- bd_t_t[-which(is.nan(cv))]
    cv <- cv[-which(is.nan(cv))]
  }

  k_hat = which(cv == min(cv))
  return(list(bd_t_t[k_hat],k_hat))
}

cv_h1_tv_cpp <- function(M = 0.1,X,FF,kmax = 100,netgd = c(0.01,1),trimp = 0.2,type = c("knn","net")){
  t = length(FF)
  MM = ceiling(M*t)-1
  w  =  abs(FF[1:(t-1)] - FF[t])
  if(length(which(duplicated(w) == T))>0){
    w  = w[-which(duplicated(w) == T)]
  }
  w = sort(w)
  wl = quantile(w,trimp)
  wu = quantile(w,1-trimp)
  ww = seq(from = wl,to = wu,length.out  = kmax+1)
  if(type == "net"){
    ww = seq(from = netgd[1],to = netgd[2],length.out  = kmax+1)
  }
  cv = cv_h1_tv_rcpp(as.matrix(X),as.matrix(FF),kmax,MM,as.matrix(ww))
  if(T %in% is.nan(cv[,2])){
    cv <- cv[-which(is.nan(cv[,2])),]
  }
  k_hat = which(cv[,2] == min(cv[,2]))
  bd_best = cv[k_hat,1]
  return(list(bd_best,k_hat,cv))
}

cv_h1_tv_fix <- function(M = 0.1,X,FF,kmax){
  cv = vector()
  bd_t_t = vector()
  t = length(FF)
  MM = ceiling(M*t)-1
  for (k in 1:kmax){
    sq_total = 0
    bd = k/(20*kmax)
    for (tt in (t-MM):t){
      a = npest_tv(FF[tt-1],bd,X[1:(tt-1),],FF[1:(tt-1)])
      cof = a
      sq = fnorm(X[tt,]- cof%*%c(1,FF[tt]))
      sq_total = sq_total + sq
    }
    bd_t_t <- c(bd_t_t,bd)
    cv <- c(cv,sq_total)
    # if(sq_total > tail(cv,1)){
    #   break
    # }
  }
  if(T %in% is.nan(cv)){
    bd_t_t <- bd_t_t[-which(is.nan(cv))]
    cv <- cv[-which(is.nan(cv))]
  }

  k_hat = which(cv == min(cv))
  return(list(bd_t_t[k_hat],k_hat))
}

K<- function(x){1/sqrt(2*pi)*exp(-x^2/2)}

mk_kernel = function(u,h2,kk,X){
  FF = X
  t = length(FF)
  mk =  sum(FF[2:t]^kk*K((FF[1:(t-1)]-u)/h2)/h2)/sum(K((FF[1:(t-1)]-u)/h2)/h2)
  return(mk)
}

cv_h2_tv <- function(M = 0.2,FF,kmax,kk){
  cv = vector()
  bd_t_t = vector()
  t = length(FF)
  MM = ceiling(M*t)-1
  for (k in 10:kmax){
    sq_total = vector()
    w  =  abs(FF[1:(t-1)] - FF[t])
    bd = max(w[which(rank(w) %in% 1:k)])
    # bd = k/100
    for (tt in (t-MM):t){
      X = FF[1:(tt-1)]
      sq = (FF[tt]^kk - mk_kernel(tail(X,1),bd,kk = kk,X))^2
      sq_total = c(sq_total,sq)
    }
    sq_total[is.nan(sq_total)] <- max(sq_total[!is.nan(sq_total)])
    sq_total = sum(sq_total)
    bd_t_t <- c(bd_t_t,bd)
    cv <- c(cv,sq_total)
  }
  if(T %in% is.nan(cv)){
    bd_t_t <- bd_t_t[-which(is.nan(cv))]
    cv <- cv[-which(is.nan(cv))]
  }
  k_hat = which(cv == min(cv))
  return(list(bd_t_t[k_hat],k_hat))
}

cv_h2_tv_cpp <- function(M = 0.2,FF,kmax,kk,trimp = 0.1){
  t = length(FF)
  MM = ceiling(M*t)-1
  w  =  abs(FF[1:(t-1)] - FF[t])
  if(length(which(duplicated(w) == T))>0){
    w  = w[-which(duplicated(w) == T)]
  }
  w = sort(w)
  wl = quantile(w,trimp)
  wu = quantile(w,1-trimp)
  ww = seq(from = wl,to = wu,length.out  = kmax)

  cv = cv_h2_tv_rcpp(MM,as.matrix(FF),kmax,kk,t,as.matrix(ww))
  if(T %in% is.nan(cv)){
    cv <- cv[-which(is.nan(cv))]
    ww <- ww[-which(is.nan(cv))]
  }
  k_hat = which(cv == min(cv))
  bd_best = ww[k_hat]
  return(list(bd_best,k_hat))
}

 # xlt = rnorm(501)
 # xnlt = matrix(rnorm(501*2),501,2)
 # et = rnorm(501)
 # zt = rnorm(501)
 # mu = exp(-zt[1:500]^2)
 # beta1 = exp(zt[1:500]^2/(1+zt[1:500]^2))
 # beta2 = exp(zt[1:500]/(1+zt[1:500]^2))
 # yt = mu + xlt[2:501] + beta1*xnlt[2:501,1] + beta2*xnlt[2:501,2] + et[2:501]
 # y = yt[1:499]
 # xl = xlt[2:500]
 # z = zt[2:500]
 # xnl = xnlt[2:500,]
 # e = et[2:500]
 # beta1[500]
 # beta2[500]
 # mu[500]
 
PLSVM = function(y,xl,xnl,z){
  y = as.matrix(y)
  xl = as.matrix(xl)
  xnl = as.matrix(xnl)
  z = as.matrix(z)
  
  t = NROW(y)
  q = NCOL(xnl)
  ql = NCOL(xl)
  
  result_np = list(bw = NA,mse = NA, beta_fit = NA)
  
  data_lm = data.frame(y = y,xl = xl, xnl = xnl)
  reg = lm(y~.,data = data_lm)
    
  e_i_lm = reg$residuals
  mu_i_lm = reg$coefficients[1]
  beta_i_lm = as.matrix(reg$coefficients[-1])
  
  result_lm = list(beta = beta_i_lm,mu = mu_i_lm,e = e_i_lm)
  
  if(ql == 0 & q > 0){
    qy = y
    qxnl = xnl
   
    #np part
    np =  try(np::npscoef(txdat = qxnl[2:t,],tydat = qy[2:t,],tzdat = z[1:(t-1),],
                          betas = TRUE,residuals = T,leave.one.out = T),silent = T)
    
    if(class(np) == "try-error"){
      beta_i = beta_i_lm
      mu_i = mu_i_lm 
      e_i = e_i_lm

      result_np = list(bw = NULL,mse = NULL,beta_fit = NULL,bws = NULL)  
    }else{
      npbw = np$bw
      
      mu_np <- predict(np, exdat = t(rep(0,q)),ezdat = z[t,])
      beta_np <- predict(np, exdat = diag(1,q,q),ezdat = rep(z[t,],q)) - mu_np
      np_bw = np$bw
      np_mse = np$MSE
      np_beta = np$beta
      e_i = vector()
      e_i[2:t] <- y[2:t,] - rowSums(cbind(1,xnl[2:t,])*np$beta)
      e_i[1] = e_i_lm[1]
      
      beta_i = beta_np
      mu_i = mu_np
      
      result_np = list(bw = np_bw,mse = np_mse,beta_fit = np_beta,bws = np$bws)  
    }
    
  }
  
  if(ql > 0 & q > 0){
    #lm part 
    Q = diag(1,t,t) - xnl%*%MASS::ginv(t(xnl)%*%xnl)%*%t(xnl)
    qy = Q%*%y
    qxl = Q%*%xl
    
    beta_pl = solve(t(qxl)%*%Q%*%qxl)%*%t(qxl)%*%Q%*%qy
    
    y_star = y - xl%*%beta_pl
    
    #np part
    np =  try(np::npscoef(txdat = xnl[2:t,],tydat = y_star[2:t,],tzdat = z[1:(t-1),],
                          betas = TRUE,residuals = T,leave.one.out = T),silent = T)
    
    
    if(class(np) == "try-error"){
      beta_i = beta_i_lm
      mu_i = mu_i_lm 
      e_i = e_i_lm
      
      result_np = list(bw = NULL,mse = NULL,beta_fit = NULL,bws = NULL)  
    }else{
      npbw = np$bw
      
      mu_np <- predict(np, exdat = t(rep(0,q)),ezdat = z[t,])
      beta_np <- predict(np, exdat = diag(1,q,q),ezdat = rep(z[t,],q))- mu_np
      np_bw = np$bw
      np_mse = np$MSE
      np_beta = np$beta
      e_i = vector()
      e_i[2:t] <- y[2:t,] - rowSums(cbind(1,xnl[2:t,])*np$beta) - xl[2:t,]%*%beta_pl 
      e_i[1] <- e_i_lm[1]
      
      beta_i = c(beta_pl,beta_np)
      mu_i = mu_np
      
      result_np = list(bw = np_bw,mse = np_mse,beta_fit = np_beta,bws = np$bws)
    }
    
    
  }
  
  if(q == 0){
    beta_i = beta_i_lm
    mu_i = mu_i_lm 
    e_i = e_i_lm
  }
  
  return(list(beta = beta_i,
              mu   = mu_i,
              e    = e_i,
              result_lm = result_lm,
              result_np = result_np))
}

PLSVM_BW = function(y,xl,xnl,z,bw = NULL){
  y = as.matrix(y)
  xl = as.matrix(xl)
  xnl = as.matrix(xnl)
  z = as.matrix(z)
  
  t = NROW(y)
  q = NCOL(xnl)
  ql = NCOL(xl)
  
  result_np = list(mu = NA,beta_np = NA)
  
  data_lm = data.frame(y = y,xl = xl, xnl = xnl)
  reg = lm(y~.,data = data_lm)
  
  e_i_lm = reg$residuals
  mu_i_lm = reg$coefficients[1]
  beta_i_lm = as.matrix(reg$coefficients[-1])
  
  result_lm = list(beta = beta_i_lm,mu = mu_i_lm,e = e_i_lm)
  
  if(ql == 0 & q > 0){
    #np part
    
    np =  npest_tv(z[t,],bw,y,xnl,z)
    
    mu_i = np[1,]
    
    beta_i = np[2:(q+1),]
    
    result_np = list(mu = mu_i,beta = beta_i)
  }
  
  if(ql > 0 & q > 0){
    Q = diag(1,t,t) - xnl%*%MASS::ginv(t(xnl)%*%xnl)%*%t(xnl)
    qy = Q%*%y
    qxl = Q%*%xl
    
    beta_pl = solve(t(qxl)%*%Q%*%qxl)%*%t(qxl)%*%Q%*%qy
    
    y_star = y - xl%*%beta_pl
    #np part
    np =  npest_tv(z[t],bw,y_star,xnl,z)
    
    mu_i = np[1,]
    
    beta_i = c(beta_pl,np[2:(1+q),])
    
    result_np = list(mu = mu_i,beta = beta_i)
  }
  
  if(q == 0){
    beta_i = beta_i_lm
    mu_i = mu_i_lm 
    e_i = e_i_lm
  }
  
  return(list(beta = beta_i,
              mu   = mu_i,
              result_lm = result_lm,
              result_np = result_np))
}

# t = 500
# xlt = rnorm(t+1)
# xnlt = matrix(rnorm((t+1)*2),t+1,2)
# et = rnorm(t+1)
# zt = rnorm(t+1)
# mu = exp(-zt[1:t]^2)
# beta1 = exp(zt[1:t]^2/(1+zt[1:t]^2))
# beta2 = exp(zt[1:t]/(1+zt[1:t]^2))
# yt = mu + xlt[2:(t+1)] + beta1*xnlt[2:(t+1),1] + beta2*xnlt[2:(t+1),2] + et[2:(t+1)]
# yf = y[t]
# y = yt[1:(t-1)]
# xl = xlt[2:t]
# z = zt[2:t]
# xnl = xnlt[2:t,]
# e = et[2:t]
# beta1[t]
# beta2[t]
# mu[t]

# ff = cbind(xnl,xl)
# x = y
# z = z

# cv_model_sel(x,ff,z,0.1)

cv_model_sel = function(x,ff,z,M = 0.1,bw.sel = c("dmz","cv"),...){
  qq = NCOL(ff)
  t = NROW(x)
  
  x = as.matrix(x)
  ff = as.matrix(ff)
  z = as.matrix(z)
  
  t_t = ceiling(M*t)
  t0 = t - t_t
  
  con_tol <- list()
  
  mse_mat = matrix(99999,choose(qq,ceiling(qq/2)),qq)
  
  for (ii in 1:qq) {
  cb =  gtools::combinations(qq,ii)
  n_cb = NROW(cb)
  
  mse_tol <- vector()
  for (jj in 1:n_cb){
    cb_k = cb[jj,]
    
    if(bw.sel == "cv"){
      inl = PLSVM(x,xl = ff[,-cb_k],xnl = ff[,cb_k], z = z)
      bw = inl$result_np$bw
    }
    if(bw.sel == "dmz"){
      bw = sd(x)*1.06*t^(-0.2)
    }
    
    
    mse = vector()
    
    for (tt in 1:t_t) {
      
      x_rw = x[tt:(tt+t0-1),]
      ffl_rw = ff[tt:(tt+t0-1),-cb_k]
      ffnl_rw = ff[tt:(tt+t0-1),cb_k]
      z_rw = z[tt:(tt+t0-1),]
      
      rw = PLSVM_BW(x_rw,xl = ffl_rw,xnl = ffnl_rw, z = z_rw,bw = bw)
      
      mse[tt] = (x[tt+t0,] - sum(c(ff[tt+t0,-cb_k],ff[tt+t0,cb_k])*rw$beta) + rw$mu)^2
      
    } 

   mse_mat[jj,ii] = sum(mse)
  }  
  
  con_tol[[ii]]  <- cb

  }
  
  mse_lm = 0
  
  for (tt in 1:t_t) {
    
    x_rw = x[tt:(tt+t0-1),]
    ffl_rw = ff[tt:(tt+t0-1),cb_k]
    ffnl_rw = ff[tt:(tt+t0-1),-cb_k]
    z_rw = z[tt:(tt+t0-1),]
    
    rw = PLSVM_BW(x_rw,xl = ffl_rw,xnl = ffnl_rw, z = z_rw,bw = bw)
    
    mse_lm = mse_lm + (x[tt+t0,] - sum(c(ff[tt+t0,cb_k],ff[tt+t0,-cb_k])*rw$beta) + rw$mu)^2
    
  }
   
  # if(mse_lm < min(mse_mat)){
  #   return(list(best_model = NA, mse_mat = mse_mat,mse_lm = mse_lm)) 
  # }else{
    row = which(mse_mat == min(mse_mat),arr.ind = TRUE)[1]
    col = which(mse_mat == min(mse_mat),arr.ind = TRUE)[2]
    
    best_model = con_tol[[col]][row,]
    
    return(list(best_model = best_model, mse_mat = mse_mat,mse_lm = mse_lm))
  # }
  
 
 
}


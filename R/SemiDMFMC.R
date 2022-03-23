#' [S]emi-parametirc [D]ynamic [M]ulti-[F]actor [M]ulti-[C]umulant estimation 
#'
#' @param X A matrix or data frame with t rows (samples) and n columns (variables).
#' @param ff The factors drive the data X.
#' @param Z The variable drive factor loadings.
#' @param SemiFL Logical. If \code{TRUE}, the factor loadings will be semi-parametric function.
#' @param sel.bw Bandwith selection method.
#' @param con_residual Logical, whether the moments of the residuals are time-varying or not.
#' @param best_model A vector. It gives which factor loadings are semi-parametric.
#' @param Penalty Only used for multi-factor model. Penalty regression to estimate factor loadings.
#' @param tgc.type the time varying type of TGC distribution.
#' @param rep_sim Repeat times for estimating time-varying parameters.
#' @param ... Any other passthru parameters
#' @return Estimated covariance, co-skewness, co-exkurtosis, co-kurtosis and the results of regression.
#' 
#' @examples
# n = 5;t = 500;k = 2;ff_lt = as.matrix(rnorm(t+2));ff_nlt = as.matrix(rnorm(t+2));
# et = matrix(rnorm(n*(t+2)),t+2,n);
# mu =   exp(-1/(1+Zt[1:(t+1)]^2));
# beta = exp((Zt[1:(t+1)]^2)/(1+Zt[1:(t+1)]^2));
# X1t = mu +  beta*(ff_nlt[2:(t+2),]) + 0.5*(ff_lt[2:(t+2),]) + (et[2:(t+2),1]);
# X2t = mu + beta*(ff_nlt[2:(t+2),]) + 1*(ff_lt[2:(t+2),]) + (et[2:(t+2),2]);
# Xt = cbind(X1t,X2t)
# X = Xt[1:t,];ff_nl = ff_nlt[2:(t+1)];ff_l = ff_lt[2:(t+1)];Z = Zt[2:(t+1)];
# con = SemiDMFMC(X,ff_l,ff_nl,Z)


SemiDMFMC = function(X,ff, Z = NULL, SemiFL = F,sel.bw = "uni",
                     con_residual = T,best_model = NULL,
                     Penalty = c("NONE","LASSO","MCP","SCAD"),
                     factor.control = list(var.model = "sGARCH", var.targeting = F, var.distribution = "sged", 
                                           Corr.Struture = c("ica", "dcc", "copula"), dcc.model = "DCC", 
                                           copula.model = list(copula = "mvt", method = "ML", time.varying = FALSE, 
                                                               transformation = "spd"), tgc.type = "leverage", tgc.targeting = F, 
                                           mean.model = list(armaOrder = c(1, 1)), CTGC = FALSE, rep_sim_f = 10),
                     eps.control    = list(var.model = "sGARCH", var.targeting = F, var.distribution = "sged",tgc.type = "leverage", 
                                           tgc.targeting = F, mean.model = list(armaOrder = c(0,0)), rep_sim_e = 10),...){
  n  = NCOL(X)
  t  = NROW(X)
  
  qq = NCOL(ff)
  
  ff = as.matrix(ff)
  X = as.matrix(X)
  
  if(!is.null(Z)){
    Z = as.matrix(Z)
  }
  
  # factor loading estimation
  result_npscoef = NULL
  result_lm = NULL
  
  if(SemiFL == T){
    
    print("Estimating the [Factor Loadings] with [Semi-parametric] function")
    
    if(sel.bw == "uni"){
      
      reg = lm(as.matrix(X) ~ as.matrix(ff))
      e_i_lm = reg$residuals
      mu_i_lm = reg$coefficients[1,]
      beta_i_lm = as.matrix(reg$coefficients[-1,])
      
      if((qq) > 1){
        beta_i_lm = t(reg$coefficients[-1,])
      }
      
      # one by one bw estimation
      np_bw = vector()
      np_mse = vector()
      np_beta = list()
      e_i_np  = matrix(NA,t,n)
      
      mu_np = vector()
      beta_np = matrix(NA,n,qq)
      
      if(NCOL(Z) < NCOL(X)){
        Z <- Z%*%t(rep(1,n)) 
      }
      
      for (ggg in 1:n){
        if(is.null(best_model)){
          ff = as.matrix(ff)
          ff_l = as.matrix(ff[,-c(1:qq)])
          ff_nl = as.matrix(ff[,1:qq])
          
          q  = NCOL(ff_nl)
          ql = NCOL(ff_l)
        }else{
          ff = as.matrix(ff)
          ff_l = as.matrix(ff[,-best_model[[ggg]]])
          ff_nl = as.matrix(ff[,best_model[[ggg]]])
          
          q  = NCOL(ff_nl)
          ql = NCOL(ff_l)
        }
        
        plsvm_con = PLSVM(y = X[,ggg],xl = ff_l,xnl = ff_nl, z = Z[,ggg])
        
        npbw = plsvm_con$result_np$bw
        
        mu_np[ggg] <- plsvm_con$mu
        beta_np[ggg,] <- plsvm_con$beta
        np_bw[ggg] = plsvm_con$result_np$bw
        np_mse[ggg] = plsvm_con$result_np$mse
        np_beta[[ggg]] = plsvm_con$result_np$beta_fit
        e_i_np[,ggg] <- plsvm_con$e
      }
      
      R2_ratio = mean(e_i_np^2)/mean(e_i_lm^2)
      
      result_npscoef = list(mu = mu_np,beta = beta_np, 
                            e =  e_i_np,
                            R2_ratio = R2_ratio, np_bw = np_bw,
                            np_beta_fit = np_beta)
      
      result_lm = list(mu = mu_i_lm,beta = beta_i_lm, 
                       e  = e_i_lm)
      
    }
    
    
    
    # if(sel.bw == "multi"){ #need improve
    # 
    #   h1_est = cv_h1_tv_cpp(M = 0.1,X,ff,kmax = 100,type = "knn")
    #   h1 = h1_est[[1]]
    # 
    #   est = npest_tv_rcpp(sp500_sample_i[t],h1,data_sample_i[1:t,],as.matrix(sp500_sample_i[1:t]))
    #   mu_np <- est[,1]
    #   beta_np <- est[,2]
    # 
    #   e_i_np  = matrix(NA,t,n)
    # 
    #   est_cof = apply(as.matrix(sp500_sample_i[1:(t-1)]),1,function(x){npest_tv_rcpp(x,h1,data_sample[1:t,],as.matrix(sp500_sample_i[1:t]))[,1:2]})
    # 
    #   mu_fit = t(est_cof[1:n,])
    #   beta_fit = t(est_cof[(n+1):(2*n),])
    # 
    #   e_i_np[2:t,]  = data_sample_i[2:t,] - mu_fit - beta_fit*sp500_sample_i[2:t]
    #   e_i_np[1,] <- e_i[1,]
    # 
    # }
    
    
  }else{
    
    if(NCOL(ff) == 1){
      
      print("Estimating the [Factor Loadings] of [Factors]")
      
      reg = lm(as.matrix(X) ~ as.matrix(ff))
      e_i_lm = reg$residuals
      mu_i_lm = reg$coefficients[1,]
      beta_i_lm = reg$coefficients[-1,]
    }
    
    if(NCOL(ff) > 1){

      if(Penalty == "NONE"){
        
        print("Estimating the [Factor Loadings] of [Factors]")
        
        reg = lm(as.matrix(X) ~ as.matrix(ff))
        e_i_lm = reg$residuals
        mu_i_lm = reg$coefficients[1,]
        beta_i_lm = t(reg$coefficients[-1,])
      }
      if(Penalty == "LASSO"){
        
        print("Estimating the [Factor Loadings] of [Multi Factors] with [Lasso]")
        
        coef_lm =  apply(as.matrix(X), 2, function(x){
          cvfit <- ncvreg::cv.ncvreg(X = ff, y = x, penalty="lasso",family = "gaussian")
          return(coef(cvfit))
        })
        e_i_lm = X - cbind(1,ff)%*%coef_lm
        mu_i_lm = coef_lm[1,]
        beta_i_lm = t(coef_lm[-1,])
      }
      if(Penalty == "SCAD"){
        
        print("Estimating the [Factor Loadings] of [Multi Factors] with [SCAD]")
        
        coef_lm =  apply(as.matrix(X), 2, function(x){
          cvfit <- ncvreg::cv.ncvreg(X = ff, y = x, penalty="SCAD",family = "gaussian")
          return(coef(cvfit))
        })
        e_i_lm = X - cbind(1,ff)%*%coef_lm
        mu_i_lm = coef_lm[1,]
        beta_i_lm = t(coef_lm[-1,])
      }
      if(Penalty == "MCP"){
        
        print("Estimating the [Factor Loadings] of [Multi Factors] with [MCP]")
        
        coef_lm =  apply(as.matrix(X), 2, function(x){
          cvfit <- ncvreg::cv.ncvreg(X = ff, y = x, penalty="MCP",family = "gaussian")
          return(coef(cvfit))
        })
        e_i_lm = X - cbind(1,ff)%*%coef_lm
        mu_i_lm = coef_lm[1,]
        beta_i_lm = t(coef_lm[-1,])
      }
    }
    
    result_lm = list(mu = mu_i_lm,beta = beta_i_lm, 
                     e  = e_i_lm)
  }
  
  if(NCOL(ff) > 1){
    
    print("Estimating the higher-order [Time-Varying] structure of [Multi Factors]")
    
    var.model        =  factor.control$var.model
    var.targeting    =  factor.control$var.targeting
    var.distribution =  factor.control$var.distribution
    Corr.Struture    =  factor.control$Corr.Struture
    dcc.model        =  factor.control$dcc.model
    copula.model     =  factor.control$copula.model
    tgc.type         =  factor.control$tgc.type
    tgc.targeting    =  factor.control$tgc.targeting
    mean.model       =  list(armaOrder = c(0, 0))
    CTGC             =  factor.control$CTGC
    rep_sim_f        =  factor.control$rep_sim_f
    
    result_factors   =  MFTGC_est(ff,var.model=var.model,var.targeting=var.targeting,var.distribution=var.distribution,
                                     Corr.Struture=Corr.Struture,dcc.model=dcc.model,copula.model=copula.model,tgc.type=tgc.type,
                                     tgc.targeting=tgc.targeting,mean.model=mean.model,CTGC=CTGC,rep_sim=rep_sim_f)

    factor_moments = result_factors[[4]]
    
  }
  
  if(NCOL(ff) == 1){
    var.model        =  factor.control$var.model
    var.targeting    =  factor.control$var.targeting
    var.distribution =  factor.control$var.distribution
    tgc.type         =  factor.control$tgc.type
    tgc.targeting    =  factor.control$tgc.targeting
    mean.model       =  factor.control$mean.model
    CTGC             =  factor.control$CTGC
    rep_sim_f        =  factor.control$rep_sim_f
    
    print("Estimating the higher-order [Time-Varying] structure of [Single Factor]")
    
    est_tgc = TGC_est(ff,var.model = var.model,var.targeting = var.targeting, 
                      var.distribution = var.distribution, tgc.type = tgc.type, 
                      tgc.targeting = tgc.targeting, mean.model = mean.model,
                      CTGC = CTGC,rep_sim_f = rep_sim_f)
    
    mu_f_fore = est_tgc$result_moment$mm.fore[1]
    var_f_fore = est_tgc$result_moment$mm.fore[2]
    skew_f_fore = est_tgc$result_moment$mm.fore[3]
    kurt_f_fore = est_tgc$result_moment$mm.fore[4]
    
    factor_moments = list(var_f_fore,skew_f_fore,kurt_f_fore)
    result_factors = list(result_snp = est_tgc,factor_moments = list(var_f_fore,skew_f_fore,kurt_f_fore))
  }
  
  #residual moment estimation
  result_residuals_tv_semi = NULL
  result_residuals_tv_lm = NULL
  
  result_residuals_semi = NULL
  result_residuals_lm = NULL
  
  if(con_residual == T){
    
    print("Estimating the higher-order [Constant] structure of [Idiosyncratic Errors]")
    
    var.model        =  eps.control$var.model
    var.targeting    =  eps.control$var.targeting
    var.distribution =  eps.control$var.distribution
    tgc.type         =  eps.control$tgc.type
    tgc.targeting    =  eps.control$tgc.targeting
    mean.model       =  eps.control$mean.model
    rep_sim_e        =  eps.control$rep_sim_e

    e_var_fore_mf_con <- vector()
    e_skew_fore_mf_con <- vector()
    e_kurt_fore_mf_con <- vector()
    
    con_residuals <- list()
    
    for (kk in 1:n){
      
      e_ik = e_i_lm[,kk]
      
      est_ei_snp = TGC_est(e_ik,var.model = var.model, var.targeting = var.targeting, 
                                      var.distribution = var.distribution, tgc.type = tgc.type, 
                                      tgc.targeting = tgc.targeting, mean.model = mean.model, CTGC = TRUE, rep_sim = rep_sim_e)
      
      con_residuals[[kk]] <- est_ei_snp
      
      e_var_fore_mf_con[kk] <- est_ei_snp$result_moment$mm.fore[2]
      e_skew_fore_mf_con[kk] <- est_ei_snp$result_moment$mm.fore[3]
      e_kurt_fore_mf_con[kk] <- est_ei_snp$result_moment$mm.fore[4]
      
    }
    
    residual_moments = list(e_var_fore_mf_con,e_skew_fore_mf_con,e_kurt_fore_mf_con) 
    result_residuals_semi = NULL
    result_residuals_lm = list(residual_snp = con_residuals, 
                               residual_moments = residual_moments)
    
    
    if(SemiFL == T){
      
      print("Estimating the higher-order [Constant] structure of [Idiosyncratic Errors] for [Semi Factor Loadings]")

      e_var_fore_smf_con <- vector()
      e_skew_fore_smf_con <- vector()
      e_kurt_fore_smf_con <- vector()
      
      con_residuals_semi <- list()
      
      for (kk in 1:n){
        e_ik = e_i_np[-1,kk]
        
        est_ei_snp = TGC_est(e_ik,var.model = var.model, var.targeting = var.targeting, 
                             var.distribution = var.distribution, tgc.type = tgc.type, 
                             tgc.targeting = tgc.targeting, mean.model = mean.model, CTGC = TRUE, rep_sim = rep_sim_e)
        
        con_residuals_semi[[kk]] <- est_ei_snp
        
        e_var_fore_smf_con[kk] <- est_ei_snp$result_moment$mm.fore[2]
        e_skew_fore_smf_con[kk] <- est_ei_snp$result_moment$mm.fore[3]
        e_kurt_fore_smf_con[kk] <- est_ei_snp$result_moment$mm.fore[4]
        
      }
      
      residual_moments = list(e_var_fore_smf_con,e_skew_fore_smf_con,e_kurt_fore_smf_con)
      result_residuals_semi = list(residual_snp = con_residuals_semi, 
                                   residual_moments = residual_moments)
    }
    
  }else{
    
    print("Estimating the higher-order [Time-Varying] structure of [Idiosyncratic Errors]")
    
    var.model        =  eps.control$var.model
    var.targeting    =  eps.control$var.targeting
    var.distribution =  eps.control$var.distribution
    tgc.type         =  eps.control$tgc.type
    tgc.targeting    =  eps.control$tgc.targeting
    mean.model       =  eps.control$mean.model
    rep_sim_e        =  eps.control$rep_sim_e
    
    
    e_var_fore_mf_tv <- vector()
    e_skew_fore_mf_tv <- vector()
    e_kurt_fore_mf_tv <- vector()
    
    con_residuals_tv <- list()
    
    for (kk in 1:n){
      
      e_ik = e_i_lm[,kk]
      
      est_ei_snp = TGC_est(e_ik,var.model = var.model, var.targeting = var.targeting, 
                           var.distribution = var.distribution, tgc.type = tgc.type, 
                           tgc.targeting = tgc.targeting, mean.model = mean.model, CTGC = FALSE, rep_sim = rep_sim_e)
      
      con_residuals_tv[[kk]] <- est_ei_snp
      
      e_var_fore_mf_tv[kk] <- est_ei_snp$result_moment$mm.fore[2]
      e_skew_fore_mf_tv[kk] <- est_ei_snp$result_moment$mm.fore[3]
      e_kurt_fore_mf_tv[kk] <- est_ei_snp$result_moment$mm.fore[4]
      
    }
    
    residual_moments = list(e_var_fore_mf_tv,
                            e_skew_fore_mf_tv,
                            e_kurt_fore_mf_tv)
    
    result_residuals_tv_semi = NULL
    
    result_residuals_tv_lm = list(residual_snp = con_residuals_tv, 
                                  residual_moments = residual_moments)
    
    if(SemiFL == T){
      
      print("Estimating the higher-order [Time-Varying] structure of [Idiosyncratic Errors] for [Semi Factor Loadings]")
      
      var.model        =  eps.control$var.model
      var.targeting    =  eps.control$var.targeting
      var.distribution =  eps.control$var.distribution
      tgc.type         =  eps.control$tgc.type
      tgc.targeting    =  eps.control$tgc.targeting
      mean.model       =  eps.control$mean.model
      rep_sim_e        =  eps.control$rep_sim_e
      
      e_var_fore_smf_tv <- vector()
      e_skew_fore_smf_tv <- vector()
      e_kurt_fore_smf_tv <- vector()
      
      con_residuals_tv_semi <- list()
      
      for (kk in 1:n){
        e_ik = e_i_np[-1,kk]
        
        est_ei_snp = TGC_est(e_ik,var.model = var.model, var.targeting = var.targeting, 
                             var.distribution = var.distribution, tgc.type = tgc.type, 
                             tgc.targeting = tgc.targeting, mean.model = mean.model, CTGC = FALSE, rep_sim = rep_sim_e)
        
        con_residuals_tv_semi[[kk]] <- est_ei_snp
        
        e_var_fore_smf_tv[kk] <- est_ei_snp$result_moment$mm.fore[2]
        e_skew_fore_smf_tv[kk] <- est_ei_snp$result_moment$mm.fore[3]
        e_kurt_fore_smf_tv[kk] <- est_ei_snp$result_moment$mm.fore[4]
        
      }
      residual_moments = list(e_var_fore_smf_tv,e_skew_fore_smf_tv,e_kurt_fore_smf_tv)
      result_residuals_tv_semi = list(residual_snp = con_residuals_tv_semi, 
                                      residual_moments = residual_moments)
      
    }
    
  }
  
  #semi-mf-tvsnp
  
  result_lm_con = NULL
  result_semi_con = NULL
  result_lm_tv = NULL
  result_semi_tv = NULL
  
  if(con_residual == T){
    
    result_lm_con = list(beta = beta_i_lm, mm_factor = factor_moments, mm_eps = residual_moments)
    
    if(SemiFL == T){
      
      result_semi_con = list(beta = beta_np, mm_factor = factor_moments, mm_eps = residual_moments)  
      
    }
    
  }else{
    
    result_lm_tv = list(beta = beta_i_lm, mm_factor = factor_moments, mm_eps = residual_moments)
    
    if(SemiFL == T){
      
      result_semi_tv = list(beta = beta_np, mm_factor = factor_moments, mm_eps = residual_moments)
      
    }
    
    
  }
  
  result = list(MM_semi_tv     =   result_semi_tv,
                MM_semi_con    =   result_semi_con,
                MM_lm_tv       =   result_lm_tv,
                MM_lm_con      =   result_lm_con,
                result_lm      =   result_lm,
                result_np      =   result_npscoef,
                result_ff      =   result_factors,
                result_res     =   list(result_residuals_tv_semi,
                                        result_residuals_semi,
                                        result_residuals_tv_lm,
                                        result_residuals_lm))
  
  print("SemiDMFMC completed")
  
  return(result)
  
}

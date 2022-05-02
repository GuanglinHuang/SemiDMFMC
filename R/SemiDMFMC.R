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
#' @param factor.control The parameters of factors' moments estimation.
#' @param eps.control The parameters of errors' moments estimation.
#' @param ... Any other passthru parameters.
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


SemiDMFMC = function(X,ff, Z = NULL, SemiFL = F,FL.type = c("VCM","SIVCM"),
                     con_residual = T,best_model = NULL,
                     Penalty = c("NONE","LASSO","MCP","SCAD"),
                     factorloading.control = list(alpha.type = 1, mm = 1, K = 15,
                                       Kmin = 10, Kmax = 20, np.eps = 10^-4, n.sim = 100, n.restarts = 5, 
                                       np.itermax = 500, trim = 0.05, M.B = FALSE,LB = 0, UB = 1, trace = 1),
                     factor.control = list(var.model = "sGARCH", var.targeting = F, var.distribution = "sged", 
                                           Corr.Struture = c("ica", "dcc", "copula"), dcc.model = "DCC", 
                                           copula.model = list(copula = "mvt", method = "ML", time.varying = TRUE, 
                                                               transformation = "spd"), tgc.type = "leverage", tgc.targeting = F, 
                                           mean.model = list(armaOrder = c(1, 1)), CTGC = FALSE, rep_sim_f = 10, n.sim_f = 5000),
                     error.control    = list(var.model = "sGARCH", var.targeting = F, var.distribution = "sged",tgc.type = "leverage", 
                                           tgc.targeting = F, mean.model = list(armaOrder = c(0,0)), rep_sim_e = 5, n.sim_e = 1000),...){
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
    
    cat("        ######################################################################################","\n",
        "       ########## Estimating the [Factor Loadings] with [Semi-parametric] function ##########","\n",
        "       ######################################################################################","\n",
        "\n")
    
    if(FL.type == "VCM"){
      
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
      
      e_i_np     = result_npscoef$e
      beta_np    = result_npscoef$beta
      u_np       = result_npscoef$mu
      alpha_np   = NA
      
    }
    
    
    
    if(FL.type == "SIVCM"){

      reg = lm(as.matrix(X) ~ as.matrix(ff))
      e_i_lm = reg$residuals
      mu_i_lm = reg$coefficients[1,]
      beta_i_lm = as.matrix(reg$coefficients[-1,])
      
      if((qq) > 1){
        beta_i_lm = t(reg$coefficients[-1,])
      }
      
      if(qq == 1){
        
        stop("Single Index Varying Coefficient Model is only used for Multi-factor model!!")
      
      }
 
      #control list
      alpha.type = factorloading.control$alpha.type
      mm         = factorloading.control$mm
      KKK        = factorloading.control$K
      Kmin       = factorloading.control$Kmin
      Kmax       = factorloading.control$Kmax
      np.eps     = factorloading.control$np.eps
      np.itermax = factorloading.control$np.itermax
      trim       = factorloading.control$trim
      M.B        = factorloading.control$M.B
      n.sim      = factorloading.control$n.sim
      n.restarts = factorloading.control$n.restarts
      LB.FL      = factorloading.control$LB
      UB.FL      = factorloading.control$UB
      trace.FL   = factorloading.control$trace
      if(alpha.type == 1){
        
        lmInitial = lm(as.matrix(X) ~ as.matrix(Z))
        
        alphaInitial = rescale(rowMeans(lmInitial$coefficients)[-1])

        alpha_inl = alphaInitial
        
        result_np = coef.bspline(yy = X,xx = ff, zz = Z, alpha_inl = alpha_inl, mm = mm, K = KKK,
                                 Kmin = Kmin, Kmax = Kmax, eps = np.eps, itermax = np.itermax, 
                                 n.sim = n.sim, n.restarts = n.restarts, trim = trim, Maxisbest = M.B,
                                 LB = LB.FL, UB = UB.FL, trace = trace.FL)
        
      }
      
      if(alpha.type == 0){
        alpha_inl = NULL
        
        result_np = coef.bspline(yy = X,xx = ff, zz = Z, alpha_inl = alpha_inl, mm = mm, K = KKK,
                                 Kmin = Kmin, Kmax = Kmax, eps = np.eps, itermax = np.itermax, 
                                 n.sim = n.sim, n.restarts = n.restarts, trim = trim, Maxisbest = M.B,
                                 LB = LB.FL, UB = UB.FL)
      }
      
      if(alpha.type == 2){

        result_np = beta.bspline(yy = X,xx = ff, zz = Z, alpha = alpha0, mm = mm, K = KKK,
                                 Kmin = Kmin, Kmax = Kmax, trim = trim, Maxisbest = M.B)
        
        
      }
      

      result_npscoef = result_np
      
      result_lm = list(mu = mu_i_lm,beta = beta_i_lm, 
                       e  = e_i_lm)
      
      e_i_np     = result_npscoef$eps_est
      beta_np    = result_npscoef$beta_est
      u_np       = result_npscoef$mu_est
      alpha_np   = result_npscoef$alpha_est
      
      
    }
    
    
  }else{
    
    if(NCOL(ff) == 1){
      
      
      cat("        ######################################################################################","\n",
          "       ##########       Estimating the [Factor Loadings] of [Single Factor]        ##########","\n",
          "       ######################################################################################","\n",
          "\n")
      
      reg = lm(as.matrix(X) ~ as.matrix(ff))
      e_i_lm = reg$residuals
      mu_i_lm = reg$coefficients[1,]
      beta_i_lm = reg$coefficients[-1,]
    }
    
    if(NCOL(ff) > 1){

      if(Penalty == "NONE"){
        
        cat("        ######################################################################################","\n",
            "       ##########       Estimating the [Factor Loadings] of [Multi-Factors]        ##########","\n",
            "       ######################################################################################","\n",
            "\n")
        
        reg = lm(as.matrix(X) ~ as.matrix(ff))
        e_i_lm = reg$residuals
        mu_i_lm = reg$coefficients[1,]
        beta_i_lm = t(reg$coefficients[-1,])
      }
      if(Penalty == "LASSO"){

        cat("        ######################################################################################","\n",
            "       ########## Estimating the [Factor Loadings] of [Multi Factors] with [Lasso] ##########","\n",
            "       ######################################################################################","\n",
            "\n")

        coef_lm =  apply(as.matrix(X), 2, function(x){
          cvfit <- ncvreg::cv.ncvreg(X = ff, y = x, penalty="lasso",family = "gaussian")
          return(coef(cvfit))
        })
        e_i_lm = X - cbind(1,ff)%*%coef_lm
        mu_i_lm = coef_lm[1,]
        beta_i_lm = t(coef_lm[-1,])
      }
      if(Penalty == "SCAD"){
        
        cat("        ######################################################################################","\n",
            "       ########## Estimating the [Factor Loadings] of [Multi Factors] with [SCAD]  ##########","\n",
            "       ######################################################################################","\n",
            "\n")
        
        coef_lm =  apply(as.matrix(X), 2, function(x){
          cvfit <- ncvreg::cv.ncvreg(X = ff, y = x, penalty="SCAD",family = "gaussian")
          return(coef(cvfit))
        })
        e_i_lm = X - cbind(1,ff)%*%coef_lm
        mu_i_lm = coef_lm[1,]
        beta_i_lm = t(coef_lm[-1,])
      }
      if(Penalty == "MCP"){
        
        cat("        ######################################################################################","\n",
            "       ##########  Estimating the [Factor Loadings] of [Multi Factors] with [MCP]  ##########","\n",
            "       ######################################################################################","\n",
            "\n")
        
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
    
    cat("        ######################################################################################","\n",
        "       ##########    Estimating the [Time-Varying] structure of [Multi-Factors]    ##########","\n",
        "       ######################################################################################","\n",
        "\n")

    
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
    n.sim_f          =  factor.control$n.sim_f
      
    result_factors   =  MFTGC_est(ff,var.model=var.model,var.targeting=var.targeting,var.distribution=var.distribution,
                                     Corr.Struture=Corr.Struture,dcc.model=dcc.model,copula.model=copula.model,tgc.type=tgc.type,
                                     tgc.targeting=tgc.targeting,mean.model=mean.model,CTGC=CTGC,rep_sim=rep_sim_f,n.sim = n.sim_f)

    factor_moments = result_factors[[4]]
    
  }
  
  if(NCOL(ff) == 1){
    var.model        =  factor.control$var.model
    var.targeting    =  factor.control$var.targeting
    var.distribution =  factor.control$var.distribution
    tgc.type         =  factor.control$tgc.type
    tgc.targeting    =  factor.control$tgc.targeting
    mean.model       =  factor.control$mean.model
    rep_sim_f        =  factor.control$rep_sim_f
    n.sim_f          =  factor.control$n.sim_f
    cat("        ######################################################################################","\n",
        "       ##########    Estimating the [Time-Varying] structure of [Single-Factors]   ##########","\n",
        "       ######################################################################################","\n",
        "\n")
    
    est_tgc = TGC_est(ff,var.model = var.model,var.targeting = var.targeting, 
                      var.distribution = var.distribution, tgc.type = tgc.type, 
                      tgc.targeting = tgc.targeting, mean.model = mean.model,
                      CTGC = T,rep_sim = rep_sim_f,n.sim = n.sim_f)
    
    mu_f_fore = est_tgc$result_moment$mm.fore[1]
    var_f_fore = est_tgc$result_moment$mm.fore[2]
    skew_f_fore = est_tgc$result_moment$mm.fore[3]
    kurt_f_fore = est_tgc$result_moment$mm.fore[4]
    
    factor_moments = list(var_f_fore,skew_f_fore,kurt_f_fore)
    result_factors = est_tgc
  }
  
  #residual moment estimation
  result_residuals_tv_semi = NULL
  result_residuals_tv_lm = NULL
  
  result_residuals_semi = NULL
  result_residuals_lm = NULL
  
  if(con_residual == T){
    
    cat("        ######################################################################################","\n",
        "       ##########  Estimating the [Constant] structure of [Idiosyncratic Errors]   ##########","\n",
        "       ######################################################################################","\n",
        "\n")

    var.model        =  error.control$var.model
    var.targeting    =  error.control$var.targeting
    var.distribution =  error.control$var.distribution
    tgc.type         =  error.control$tgc.type
    tgc.targeting    =  error.control$tgc.targeting
    mean.model       =  error.control$mean.model
    rep_sim_e        =  error.control$rep_sim_e
    n.sim_e          =  error.control$n.sim_e
    
    e_var_fore_mf_con <- vector()
    e_skew_fore_mf_con <- vector()
    e_kurt_fore_mf_con <- vector()
    
    con_residuals <- list()
    
    for (kk in 1:n){
      
      e_ik = e_i_lm[,kk]
      
      est_ei_snp = TGC_est(e_ik,var.model = var.model, var.targeting = var.targeting, 
                                      var.distribution = var.distribution, tgc.type = tgc.type, 
                                      tgc.targeting = tgc.targeting, mean.model = mean.model, CTGC = T, rep_sim = rep_sim_e,n.sim = n.sim_e)
      
      con_residuals[[kk]] <- est_ei_snp
      
      e_var_fore_mf_con[kk] <- est_ei_snp$result_moment$mm.fore[2]
      e_skew_fore_mf_con[kk] <- est_ei_snp$result_moment$mm.fore[3]
      e_kurt_fore_mf_con[kk] <- est_ei_snp$result_moment$mm.fore[4]
      
    }
    
    residual_moments = list(e_var_fore_mf_con,e_skew_fore_mf_con,e_kurt_fore_mf_con) 
    result_residuals_semi = NULL
    result_residuals_lm = con_residuals
    
    
    if(SemiFL == T){
      
      cat("        ######################################################################################","\n",
          "       #####  Estimating the [Constant] structure of [Idiosyncratic Errors] for [SFL]  ######","\n",
          "       ######################################################################################","\n",
          "\n")

      e_var_fore_smf_con <- vector()
      e_skew_fore_smf_con <- vector()
      e_kurt_fore_smf_con <- vector()
      
      con_residuals_semi <- list()
      
      for (kk in 1:n){
        
        if(FL.type == "SIVCM"){
          e_ik = e_i_np[-c(1:mm),kk]
        }else{
          e_ik = e_i_np[-1,kk]
        }

        est_ei_snp = TGC_est(e_ik,var.model = var.model, var.targeting = var.targeting, 
                             var.distribution = var.distribution, tgc.type = tgc.type, 
                             tgc.targeting = tgc.targeting, mean.model = mean.model, CTGC = T, rep_sim = rep_sim_e,n.sim = n.sim_e)
        
        con_residuals_semi[[kk]] <- est_ei_snp
        
        e_var_fore_smf_con[kk] <- est_ei_snp$result_moment$mm.fore[2]
        e_skew_fore_smf_con[kk] <- est_ei_snp$result_moment$mm.fore[3]
        e_kurt_fore_smf_con[kk] <- est_ei_snp$result_moment$mm.fore[4]
        
      }
      
      residual_moments = list(e_var_fore_smf_con,e_skew_fore_smf_con,e_kurt_fore_smf_con)
      result_residuals_semi = con_residuals_semi
    }
    
  }else{
    
    cat("        ######################################################################################","\n",
        "       ######    Estimating the [Time-Varying] structure of [Idiosyncratic Errors]     ######","\n",
        "       ######################################################################################","\n",
        "\n")

    var.model        =  error.control$var.model
    var.targeting    =  error.control$var.targeting
    var.distribution =  error.control$var.distribution
    tgc.type         =  error.control$tgc.type
    tgc.targeting    =  error.control$tgc.targeting
    mean.model       =  error.control$mean.model
    rep_sim_e        =  error.control$rep_sim_e
    n.sim_e          =  error.control$n.sim_e
    
    e_var_fore_mf_tv <- vector()
    e_skew_fore_mf_tv <- vector()
    e_kurt_fore_mf_tv <- vector()
    
    con_residuals_tv <- list()
    
    for (kk in 1:n){
      
      e_ik = e_i_lm[,kk]
      
      est_ei_snp = TGC_est(e_ik,var.model = var.model, var.targeting = var.targeting, 
                           var.distribution = var.distribution, tgc.type = tgc.type, 
                           tgc.targeting = tgc.targeting, mean.model = mean.model, CTGC = FALSE, rep_sim = rep_sim_e,n.sim = n.sim_e)
      
      con_residuals_tv[[kk]] <- est_ei_snp
      
      e_var_fore_mf_tv[kk] <- est_ei_snp$result_moment$mm.fore[2]
      e_skew_fore_mf_tv[kk] <- est_ei_snp$result_moment$mm.fore[3]
      e_kurt_fore_mf_tv[kk] <- est_ei_snp$result_moment$mm.fore[4]
      
    }
    
    residual_moments = list(e_var_fore_mf_tv,
                            e_skew_fore_mf_tv,
                            e_kurt_fore_mf_tv)
    
    result_residuals_tv_semi = NULL
    
    result_residuals_tv_lm = con_residuals_tv
    
    if(SemiFL == T){
      
      cat("        ######################################################################################","\n",
          "       ###  Estimating the [Time-Varying] structure of [Idiosyncratic Errors] for [SFL]   ###","\n",
          "       ######################################################################################","\n",
          "\n")
      
      var.model        =  error.control$var.model
      var.targeting    =  error.control$var.targeting
      var.distribution =  error.control$var.distribution
      tgc.type         =  error.control$tgc.type
      tgc.targeting    =  error.control$tgc.targeting
      mean.model       =  error.control$mean.model
      rep_sim_e        =  error.control$rep_sim_e
      n.sim_e          =  error.control$n.sim_e
      
      e_var_fore_smf_tv <- vector()
      e_skew_fore_smf_tv <- vector()
      e_kurt_fore_smf_tv <- vector()
      
      con_residuals_tv_semi <- list()
      
      for (kk in 1:n){
        
        if(FL.type == "SIVCM"){
          e_ik = e_i_np[-c(1:mm),kk]
        }else{
          e_ik = e_i_np[-1,kk]
        }
        
        
        est_ei_snp = TGC_est(e_ik,var.model = var.model, var.targeting = var.targeting, 
                             var.distribution = var.distribution, tgc.type = tgc.type, 
                             tgc.targeting = tgc.targeting, mean.model = mean.model, CTGC = FALSE, rep_sim = rep_sim_e,n.sim = n.sim_e)
        
        con_residuals_tv_semi[[kk]] <- est_ei_snp
        
        e_var_fore_smf_tv[kk] <- est_ei_snp$result_moment$mm.fore[2]
        e_skew_fore_smf_tv[kk] <- est_ei_snp$result_moment$mm.fore[3]
        e_kurt_fore_smf_tv[kk] <- est_ei_snp$result_moment$mm.fore[4]
        
      }
      
      residual_moments = list(e_var_fore_smf_tv,e_skew_fore_smf_tv,e_kurt_fore_smf_tv)
      
      result_residuals_tv_semi = con_residuals_tv_semi
      
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
      
      result_semi_con = list(beta = beta_np,u = u_np, alpha = alpha_np, mm_factor = factor_moments, mm_eps = residual_moments)  
      
    }
    
  }else{
    
    result_lm_tv = list(beta = beta_i_lm, mm_factor = factor_moments, mm_eps = residual_moments)
    
    if(SemiFL == T){
      
      result_semi_tv = list(beta = beta_np,u = u_np, alpha = alpha_np, mm_factor = factor_moments, mm_eps = residual_moments)
      
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
  
  cat("        ######################################################################################","\n",
      "       #################                 SemiDMFMC Completed                #################","\n",
      "       ######################################################################################","\n",
      "\n")
  
  return(result)
  
}

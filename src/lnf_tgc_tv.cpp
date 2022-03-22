#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
double lnf_tgc_tv(SEXP ZZZ, SEXP TTT, SEXP T0){ 
  
  NumericVector ZZ(ZZZ);
  NumericVector THETA(TTT);
  
  double t10 = THETA[0];
  double t11 = THETA[1];
  double t121 = THETA[2];
  double t122 = THETA[3];
  double t13 = THETA[4];
  double t20 = THETA[5];
  double t21 = THETA[6];
  double t221 = THETA[7];
  double t222 = THETA[8];
  double t23 = THETA[9];
  
  int TT = as<int>(T0);
  double theta1_t1 = 0.0;
  double theta2_t1 = 0.0;
  double theta1_t;
  double theta2_t;
  
  double gam1;
  double gam2;
  double lam;
  
  double Ex;double Ex2;double b_theta_t;double a_theta_t;double x_t;double H3;
  double H4;double phi;
  double fff = 0.0;
  
  for(int tt = 1; tt < TT; tt++) {
     
     double zt1 = ZZ[tt-1];
     double zt = ZZ[tt];
  
     if(zt1 > 0){
       theta1_t = t10 + t11*theta1_t1 + t121*(1 + t13*abs(zt1))*zt1; 
       theta2_t = t20 + t21*theta2_t1 + t221*(1 + t23*abs(zt1))*zt1;
     }
     if(zt1 <= 0){
       theta1_t = t10 + t11*theta1_t1 + t122*(1 + t13*abs(zt1))*zt1; 
       theta2_t = t20 + t21*theta2_t1 + t222*(1 + t23*abs(zt1))*zt1;
     }

     gam1 = theta1_t/sqrt(6);
     gam2 = theta2_t/sqrt(24);
     lam = 1/(1+gam1*gam1+gam2*gam2);

     Ex = 4*lam*gam1*gam2;
     Ex2 = 6*lam*gam1*gam1 + 8*lam*gam2*gam2 + 1;

     b_theta_t = 1/sqrt(Ex2 - Ex*Ex);
     a_theta_t = - b_theta_t*Ex;
     x_t = (zt - a_theta_t)/b_theta_t;

     H3 = (x_t*x_t*x_t-3*x_t)/sqrt(6);
     H4 = (x_t*x_t*x_t*x_t-6*x_t*x_t+3)/sqrt(24);
     phi = 1 + gam1*H3 + gam2*H4;

     fff += -log(b_theta_t) + log(lam) - 0.5*x_t*x_t +  2*log(abs(phi));
     theta1_t1 = theta1_t;
     theta2_t1 = theta2_t;
   }
  return fff;
}
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
// [[Rcpp::export]]
NumericVector lnf_tgc_tv_lt(SEXP ZZZ, SEXP TTT, SEXP T0){ 
  
  NumericVector ZZ(ZZZ);
  NumericVector THETA(TTT);
  
  double t10 = THETA[0];
  double t11 = THETA[1];
  double t121 = THETA[2];
  double t122 = THETA[3];
  double t13 = THETA[4];
  double t20 = THETA[5];
  double t21 = THETA[6];
  double t221 = THETA[7];
  double t222 = THETA[8];
  double t23 = THETA[9];
  
  int TT = as<int>(T0);
  double theta1_t1 = 0.0;
  double theta2_t1 = 0.0;
  double theta1_t;
  double theta2_t;
  
  double gam1;
  double gam2;
  double lam;
  
  double Ex;double Ex2;double b_theta_t;double a_theta_t;double x_t;double H3;
  double H4;double phi;
  NumericVector fff(ZZ.length()-1);
  
  for(int tt = 1; tt < TT; tt++) {
    
    double zt1 = ZZ[tt-1];
    double zt = ZZ[tt];
    
    if(zt1 > 0){
      theta1_t = t10 + t11*theta1_t1 + t121*(1 + t13*abs(zt1))*zt1; 
      theta2_t = t20 + t21*theta2_t1 + t221*(1 + t23*abs(zt1))*zt1;
    }
    if(zt1 <= 0){
      theta1_t = t10 + t11*theta1_t1 + t122*(1 + t13*abs(zt1))*zt1; 
      theta2_t = t20 + t21*theta2_t1 + t222*(1 + t23*abs(zt1))*zt1;
    }
    
    gam1 = theta1_t/sqrt(6);
    gam2 = theta2_t/sqrt(24);
    lam = 1/(1+gam1*gam1+gam2*gam2);
    
    Ex = 4*lam*gam1*gam2;
    Ex2 = 6*lam*gam1*gam1 + 8*lam*gam2*gam2 + 1;
    
    b_theta_t = 1/sqrt(Ex2 - Ex*Ex);
    a_theta_t = - b_theta_t*Ex;
    x_t = (zt - a_theta_t)/b_theta_t;
    
    H3 = (x_t*x_t*x_t-3*x_t)/sqrt(6);
    H4 = (x_t*x_t*x_t*x_t-6*x_t*x_t+3)/sqrt(24);
    phi = 1 + gam1*H3 + gam2*H4;
    
    fff[tt-1] = -log(b_theta_t) + log(lam) - 0.5*x_t*x_t +  2*log(abs(phi));
    theta1_t1 = theta1_t;
    theta2_t1 = theta2_t;
  }
  return fff;
}

// [[Rcpp::export]]
double lnf_tgc_tv_tar(SEXP ZZZ, SEXP TTT,SEXP TTT0, SEXP T0){ 
  using namespace Rcpp ;
  
  NumericVector ZZ(ZZZ);
  NumericVector THETA(TTT);
  NumericVector THETA0(TTT0);
  
  double t10 = THETA0[0];
  double t20 = THETA0[1];
  
  double t11 = THETA[0];
  double t121 = THETA[1];
  double t122 = THETA[2];
  double t13 = THETA[3];
  
  double t21 = THETA[4];
  double t221 = THETA[5];
  double t222 = THETA[6];
  double t23 = THETA[7];
  
  int TT = as<int>(T0);
  double theta1_t1 = 0.0;
  double theta2_t1 = 0.0;
  double theta1_t;
  double theta2_t;
  
  double gam1;
  double gam2;
  double lam;
  
  double Ex;double Ex2;double b_theta_t;double a_theta_t;double x_t;double H3;
  double H4;double phi;
  double fff = 0.0;
  
  for(int tt = 1; tt < TT; tt++) {
    
    double zt1 = ZZ[tt-1];
    double zt = ZZ[tt];
    
    if(zt1 > 0){
      theta1_t = t10 + t11*theta1_t1 + t121*(1 + t13*abs(zt1))*zt1; 
      theta2_t = t20 + t21*theta2_t1 + t221*(1 + t23*abs(zt1))*zt1;
    }
    if(zt1 <= 0){
      theta1_t = t10 + t11*theta1_t1 + t122*(1 + t13*abs(zt1))*zt1; 
      theta2_t = t20 + t21*theta2_t1 + t222*(1 + t23*abs(zt1))*zt1;
    }
    
    gam1 = theta1_t/sqrt(6);
    gam2 = theta2_t/sqrt(24);
    lam = 1/(1+gam1*gam1+gam2*gam2);
    
    Ex = 4*lam*gam1*gam2;
    Ex2 = 6*lam*gam1*gam1 + 8*lam*gam2*gam2 + 1;
    
    b_theta_t = 1/sqrt(Ex2 - Ex*Ex);
    a_theta_t = - b_theta_t*Ex;
    x_t = (zt - a_theta_t)/b_theta_t;
    
    H3 = (x_t*x_t*x_t-3*x_t)/sqrt(6);
    H4 = (x_t*x_t*x_t*x_t-6*x_t*x_t+3)/sqrt(24);
    phi = 1 + gam1*H3 + gam2*H4;
    
    fff += -log(b_theta_t) + log(lam) - 0.5*x_t*x_t +  2*log(abs(phi));
    theta1_t1 = theta1_t;
    theta2_t1 = theta2_t;
  }
  return fff;
}

// [[Rcpp::export]]
NumericVector lnf_tgc_tv_tar_lt(SEXP ZZZ, SEXP TTT,SEXP TTT0, SEXP T0){ 
  using namespace Rcpp ;
  
  NumericVector ZZ(ZZZ);
  NumericVector THETA(TTT);
  NumericVector THETA0(TTT0);
  
  double t10 = THETA0[0];
  double t20 = THETA0[1];
  
  double t11 = THETA[0];
  double t121 = THETA[1];
  double t122 = THETA[2];
  double t13 = THETA[3];
  
  double t21 = THETA[4];
  double t221 = THETA[5];
  double t222 = THETA[6];
  double t23 = THETA[7];
  
  int TT = as<int>(T0);
  double theta1_t1 = 0.0;
  double theta2_t1 = 0.0;
  double theta1_t;
  double theta2_t;
  
  double gam1;
  double gam2;
  double lam;
  
  double Ex;double Ex2;double b_theta_t;double a_theta_t;double x_t;double H3;
  double H4;double phi;
  NumericVector fff(ZZ.length()-1);
  
  for(int tt = 1; tt < TT; tt++) {
    
    double zt1 = ZZ[tt-1];
    double zt = ZZ[tt];
    
    if(zt1 > 0){
      theta1_t = t10 + t11*theta1_t1 + t121*(1 + t13*abs(zt1))*zt1; 
      theta2_t = t20 + t21*theta2_t1 + t221*(1 + t23*abs(zt1))*zt1;
    }
    if(zt1 <= 0){
      theta1_t = t10 + t11*theta1_t1 + t122*(1 + t13*abs(zt1))*zt1; 
      theta2_t = t20 + t21*theta2_t1 + t222*(1 + t23*abs(zt1))*zt1;
    }
    
    gam1 = theta1_t/sqrt(6);
    gam2 = theta2_t/sqrt(24);
    lam = 1/(1+gam1*gam1+gam2*gam2);
    
    Ex = 4*lam*gam1*gam2;
    Ex2 = 6*lam*gam1*gam1 + 8*lam*gam2*gam2 + 1;
    
    b_theta_t = 1/sqrt(Ex2 - Ex*Ex);
    a_theta_t = - b_theta_t*Ex;
    x_t = (zt - a_theta_t)/b_theta_t;
    
    H3 = (x_t*x_t*x_t-3*x_t)/sqrt(6);
    H4 = (x_t*x_t*x_t*x_t-6*x_t*x_t+3)/sqrt(24);
    phi = 1 + gam1*H3 + gam2*H4;
    
    fff[tt-1] = -log(b_theta_t) + log(lam) - 0.5*x_t*x_t +  2*log(abs(phi));
    theta1_t1 = theta1_t;
    theta2_t1 = theta2_t;
  }
  return fff;
}

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector cv_h2_tv_rcpp(SEXP MM, SEXP FFF, SEXP KMAX, double KK,SEXP TT, SEXP WW){ 
  using namespace Rcpp;
  NumericVector F(FFF);
  int kmax = as<int>(KMAX);
  int mm = as<int>(MM);
  int t = as<int>(TT);
  NumericVector ww(WW);
  NumericVector cv(kmax);
  NumericVector bd_t_t(kmax);
  
  double pi = 3.14159265358979323846;
  double exp = 2.71828182845904523536;
 
  double bd;
  double kk = KK;
for(int k = 0; k < kmax; k++) {
      bd = ww[k];
      double mse = 0.0;
      
for(int tt = (t-mm-1); tt < t; tt++){
    double fz = 0.0;
    double fm = 0.0;
    double mk = 0.0;
    double Kern = 0.0;
    for(int jj = 1; jj < tt-1; jj++){
    Kern = 1/sqrt(2*pi)*pow(exp,-((F[jj-1]-F[tt-2])/bd)*((F[jj-1]-F[tt-2])/bd)/2);
    fz += (pow(F[jj],kk))*Kern/bd;
    fm += Kern/bd;
    }
    if(fm == 0){fm = 0.0000001;}
    mk = fz/fm;
    mse += (F[tt-1] - mk)*(F[tt-1] - mk);

}
     bd_t_t[k] = bd;
     cv[k] = mse;
}
   
  return cv;
}






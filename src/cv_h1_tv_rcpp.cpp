#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat cv_h1_tv_rcpp(SEXP YY, SEXP FF, SEXP KMAX, SEXP MM, SEXP BD){

  int kmax = as<int>(KMAX);
  int mm = as<int>(MM);

  arma::mat Y = as<arma::mat>(YY);
  arma::mat X = as<arma::mat>(FF);
  int t = X.n_rows;
  arma::mat Bandwidth = as<arma::mat>(BD);
  arma::mat bd_t_t(kmax,2);

for(int k = 0; k < kmax; k++) {
      double sq_total = 0.0;
      double bd = Bandwidth(k,0);
for(int tt = (t-mm-1); tt < t; tt++){

  double h = bd;
  double u = X(tt-2,0);
  double pi = 3.14159265358979323846;
  arma::mat YYY = Y.rows(0,tt-2);
  arma::mat XXX = X.rows(0,tt-2);

  int tj = XXX.n_rows;
  arma::mat YMY = YYY.rows(1,tj-1);
  arma::mat W(tj-1,tj-1);
  arma::mat XX(tj-1,4);
  arma::mat uu(tj-1,1);
  arma::mat ee(tj-1,1);
  arma::mat hh(tj-1,1);
  arma::mat two(tj-1,1);
  arma::mat pipi(tj-1,1);

  W.fill(0);uu.fill(u); ee.fill(1.0);hh.fill(h);two.fill(2.0);pipi.fill(pi);
  XX.col(0) = ee; XX.col(1) = XXX.rows(1,tj-1);
  XX.col(2) = XXX.rows(0,tj-2) - uu;
  XX.col(3) = XXX.rows(1,tj-1)%(XXX.rows(0,tj-2)- uu);
  arma::mat f = XXX.rows(0,tj-2);

  W.diag() = (ee/sqrt(two%pipi))%exp(-((uu - f)/h)%((uu - f)/h)/two)/h;

  arma::mat Beta = pinv(XX.t()*W*XX)*XX.t()*W*YMY;
  arma::mat cof = Beta.rows(0,1);
  arma::mat XXT(1,2);

  XXT(0,0) = 1.0; XXT(0,1) = X(tt-1,0);
  sq_total += norm(Y.row(tt-1) - XXT*cof,"fro");
  }
  bd_t_t(k,0) = bd;
  bd_t_t(k,1) = sq_total;
}
  return bd_t_t;
}






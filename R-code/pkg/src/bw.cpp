// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

# include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List ComputeAlphaBeta(arma::vec y, arma::vec x, arma::mat WW, arma::mat weight,
                      arma::mat Z, arma::mat B, arma::mat Bk) {
    arma::mat weightSQR = sqrt(weight);

    x.each_col() %= weightSQR;
    y.each_col() %= weightSQR;
    Z.each_col() %= weightSQR;
    WW.each_col() %= weightSQR;
    B.each_col() %= weightSQR;
    Bk.each_row() %= weightSQR.t();

    weightSQR.reset();

    arma::mat xx = x - WW * solve(WW, x);
    x.reset();
    arma::mat yy = y - WW * solve(WW, y);
    y.reset();
    arma::mat ZZ = Z - WW * solve(WW, Z);
    WW.reset();

    arma::vec Zxx = Z.t() * xx;
    arma::vec Zyy = Z.t() * yy;
    Z.reset();
    arma::colvec ZZZZ = sum(ZZ.t() % ZZ.t(), 1);

    arma::mat Alpha = (Bk * xx) / as_scalar(B.t() * xx);
    arma::mat Beta = Zyy / Zxx;
    arma::mat Gamma = Zyy / ZZZZ;
    arma::mat pi = (ZZ.t() * xx) / ZZZZ;

    return List::create(Alpha, Beta, Gamma, pi);
}
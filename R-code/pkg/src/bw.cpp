// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

# include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List ComputeAlphaBeta(arma::vec y, arma::vec x, arma::mat WW,
                      arma::mat weight, arma::mat Z, arma::mat G) {
    arma::mat weightSQR = sqrt(weight);

    x.each_col() %= weightSQR;
    y.each_col() %= weightSQR;
    Z.each_col() %= weightSQR;
    WW.each_col() %= weightSQR;
    weightSQR.reset();

    arma::vec xx = x - WW * solve(WW, x);
    x.reset();
    arma::vec yy = y - WW * solve(WW, y);
    y.reset();
    arma::mat ZZ = Z - WW * solve(WW, Z);
    WW.reset();

    arma::vec Zxx = Z.t() * xx;
    arma::vec Zyy = Z.t() * yy;
    arma::colvec ZZZZ = sum(ZZ.t() % ZZ.t(), 1);

    arma::mat Alpha = G.each_col() % Zxx;
    Alpha.each_row() /= Zxx.t() * G;

    arma::vec Beta = Zyy / Zxx;
    arma::vec Gamma = Zyy / ZZZZ;
    arma::vec pi = (ZZ.t() * xx) / ZZZZ;

    return List::create(Alpha, Beta, Gamma, pi);
}

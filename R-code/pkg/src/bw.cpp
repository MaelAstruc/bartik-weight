
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

# include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List ComputeAlphaBeta(arma::vec y, arma::vec x, arma::mat WW,
                      arma::vec weightVec, arma::mat Z, arma::mat G) {

    int n = x.n_elem;
    int nrow_G = G.n_rows;
    int ncol_G = G.n_cols;

    // TO DO: optimize memory usage
    arma::mat weight = diagmat(weightVec);

    arma::mat I;
    arma::mat M_W = I.eye(n, n) - WW * inv_sympd(WW.t() * weight * WW) * WW.t() * weight;
    arma::vec xx = M_W * x;
    arma::vec yy = M_W * y;
    arma::mat ZZ = M_W * Z;
    arma::mat Alpha(nrow_G, ncol_G);
    for (int i = 0; i < ncol_G; i++){
        arma::vec g = G.col(i);
        Alpha.col(i) = (diagmat(g) * Z.t() * weight * xx) / as_scalar(g.t() * Z.t() * weight * xx);
    }
    arma::vec Beta = (Z.t() * weight * yy) / (Z.t() * weight * xx);
    arma::vec Gamma = (Z.t() * weight * yy) / ((ZZ.t() % ZZ.t()) * weightVec);
    arma::vec pi = (ZZ.t() * weight * xx) / ((ZZ.t() % ZZ.t()) * weightVec);

    return List::create(Alpha, Beta, Gamma, pi);
}

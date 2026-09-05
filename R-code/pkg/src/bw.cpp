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

    // Compute reduced form residuals per group
    arma::mat res_0 = arma::mat(yy.n_rows, Beta.n_rows);
    for (unsigned int i = 0; i < Beta.n_rows; i++) res_0.col(i) = yy.col(0) - as_scalar(Beta.row(i)) * xx;

    // Compute first-stage residuals per group
    arma::mat res_1 = arma::mat(xx.n_rows, pi.n_rows);
    for (unsigned int i = 0; i < pi.n_rows; i++) res_1.col(i) = xx.col(0) - as_scalar(pi.row(i)) * ZZ.col(i);

    // Prepare matrices to retrieve diagonals
    arma::mat res_0_2 = res_0.t() * res_0;
    arma::mat res_1_2 = res_1.t() * res_1;
    arma::mat xzzx = Zxx * Zxx.t();
    arma::mat zz = ZZ.t() * ZZ;
    arma::mat pi_2 = pi * pi.t();

    arma::mat se_0 = sqrt(res_0_2.diag() / xzzx.diag() / (double) (xx.n_rows -  xx.n_cols));
    arma::mat se_1 = res_1_2.diag() / zz.diag() / (double) (ZZ.n_rows - ZZ.n_cols);

    arma::mat f1 = pi_2.diag() / se_1;

    return List::create(Alpha, Beta, Gamma, pi, se_0, f1, pi_2, se_1);
}

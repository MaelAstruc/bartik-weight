#' Estimating Rotemberg Weights for a Bartik "Instrument"
#'
#' `bw()` estimates the Rotemberg weights for a Bartik "instrument" outlined in
#' Goldsmith-Pinkham, Sorkin and Swift (2018) and the just-identified IV
#' estimates using local industry shares, which are the actual instruments.
#'
#' The typical Bartik "instrument" has two components: the local industry
#' shares (usually at the level of location-industry or location-year-industry),
#' and the overall industry growth rates (usually at the level of industry or
#' industry-year). The outcome variable and the endogenous variable are at the
#' level of location or location-year.
#'
#' Regarding the application in the economics of migration, the shock comes
#' from the changes in the population of migrants, the shift. The share is the
#' distribution of migrants, grouped by origin, across the locations in a
#' previous period of reference.
#'
#' Because the key variables are at three different levels, `bw()` proceeds by
#' taking three different datasets: (1) a main `master` data frame containing
#' the dependent variable (`y`), the causal variable of interest (`x`), a set of
#' control variables (`controls`), and the weighted variables (`weight`); (2) a
#' `local` data frame containing the shares (`Z`); and (3) a `global` data frame
#' containing the shift (`G`). At the moment, it is necessary to transform the
#' "local" dataset from long format to wide format.
#'
#' Given the dimensions:
#' - N the number of observations
#' - L the number of localions
#' - K the number of groups
#' - T the number of periods
#'
#' The dimensions of these datasets are :
#' - `master` : N rows
#' - `local[Z]` : N rows and K columns
#' - `global[G]` : K rows and T columns
#'
#' Be sure to match the `master` and `local` datasets rows beforehand and
#' to match the order of the `local` dataset column with the `global` dataset
#' rows.
#'
#' If the number of rows in the `local` dataset is different from N because
#' your `master` dataset is not at the localion or location-industry level, you
#' need to duplicate the `local` rows to match with the `master` rows.
#'
#' @md
#' @param master The master data frame.
#' @param y A string for outcome variable. It should be a variable in `master`.
#' @param x A string for the endogenous variable. It should be a variable in
#'   `master`.
#' @param controls A string or character vector for the control variables.
#'   `controls` are optional and should be in `master`.
#' @param weight A string for the weighted variable. `weight` is optional and
#'   should be in `master`.
#' @param local The local data frame. It should be in wide format.
#' @param Z A string or character vector for the the local industry shares.
#' @param global The global data frame.
#' @param G A string for the the overall industry growth rates.
#'
#' @return A tibble with K rows, a for each group, with :
#' \itemize{
#'   \item the `global` dataset.
#'   \item one column or multiple `alpha` for the Rotemberg weights. If $T$ is larger than 1, the columns are named `alpha_'t'` with `t` the names in @param G.
#'   \item one column `beta` for the just-identified IV estimated.
#'   \item one column `gamma` reduced form/intent-to-treat effect of the Bartik IV on Y.
#'   \item one column `pi` the first stage effect of the Bartik IV on X.
#' }
#'
#' @importFrom tibble as_tibble
#' @importFrom Rcpp sourceCpp
#' @useDynLib bartik.weight, .registration = TRUE
#'
#' @examples
#' library(bartik.weight)
#'
#' index = c("czone", "year")
#' y = "d_sh_empl_mfg"
#' x = "d_tradeusch_pw"
#' controls = c("reg_midatl", "reg_encen", "reg_wncen", "reg_satl",
#'   "reg_escen", "reg_wscen", "reg_mount", "reg_pacif", "l_sh_popedu_c",
#'   "l_sh_popfborn", "l_sh_empl_f", "l_sh_routine33", "l_task_outsource",
#'   "t2", "l_shind_manuf_cbp")
#' weight = "timepwt48"
#' Z = setdiff(names(ADH_local_wide), index)
#' G = "trade_"
#'
#' bw(ADH_master, y, x, controls, weight, ADH_local_wide, Z, ADH_global, G)
#'
#' @export
bw = function(master, y, x, controls = NULL, weight = NULL, local, Z, global, G) {
    # Checks
    stopifnot(nrow(master) == nrow(local))
    stopifnot(length(Z) == nrow(global))

    # Parsing the master file
    y = master[[y]]
    x = master[[x]]
    n = length(x)

    if (is.null(weight)) {
        weight = rep(1, n)
    } else {
        weight = master[[weight]]
    }

    if (is.null(controls)) {
        WW = matrix(1, n, 1)
    } else {
        W = as.matrix(master[controls])
        WW = cbind(W, matrix(1, n, 1))
    }

    # Parsing the local file
    Z = as.matrix(local[Z])

    # Parsing the global file
    G = as.matrix(global[G])

    # Compute the Rotemberg weights (alpha) and the just-identified coefficients (beta)
    alpha_beta = ComputeAlphaBeta(y, x, WW, weight, Z, G)

    # Prepare matrice of alphas
    alphas = alpha_beta[[1]]
    if (ncol(alphas) == 1) {
        colnames(alphas) <- "alpha"
    } else {
        colnames(alphas) <- paste0("alpha_", colnames(G))
    }

    # Return a tibble
    tibble::as_tibble(cbind(global, alphas,
                            beta = alpha_beta[[2]],
                            gamma = alpha_beta[[3]], pi = alpha_beta[[4]])
                      )
}

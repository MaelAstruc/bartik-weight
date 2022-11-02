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
#' - L the number of locations
#' - K the number of groups
#' - T the number of periods
#'
#' The dimensions of these datasets are :
#' - `master` : N rows
#' - `local[Z]` : L rows and K columns
#' - `global[G]` : K rows and T columns
#'
#' @md
#'
#' @param master The master data frame.
#' @param y A string for outcome variable. It should be a variable in `master`.
#' @param x A string for the endogenous variable. It should be a variable in
#'   `master`.
#' @param controls A string or character vector for the control variables.
#'   `controls` are optional and should be in `master`.
#' @param weight A string for the weighted variable. `weight` is optional and
#'   should be in `master`.
#' @param local The local data frame. It should be in wide format.
#' @param z The variables used in the local data frame.
#' @param global The global data frame.
#' @param g The variables used in the global data frame.
#' @param L_id The indexes used to identify and match the locations in the
#'   master and local data frames.
#' @param T_id The indexes used to identify and match the periods in the master
#'   and global data frames.
#'
#' @return A tibble with K rows, a for each group, with :
#' \itemize{
#'   \item the `global` dataset.
#'   \item one column `alpha` for the Rotemberg weights.
#'   \item one column `beta` for the just-identified IV estimated.
#'   \item one column `gamma` reduced form/intent-to-treat effect of the Bartik IV on Y.
#'   \item one column `pi` the first stage effect of the Bartik IV on X.
#' }
#'
#' @importFrom tibble as_tibble
#'
#' @examples
#' library(bartik.weight)
#'
#' index <- c("czone", "year")
#' y <- "d_sh_empl_mfg"
#' x <- "d_tradeusch_pw"
#' controls <- c("reg_midatl", "reg_encen", "reg_wncen", "reg_satl",
#'   "reg_escen", "reg_wscen", "reg_mount", "reg_pacif", "l_sh_popedu_c",
#'   "l_sh_popfborn", "l_sh_empl_f", "l_sh_routine33", "l_task_outsource",
#'   "t2", "l_shind_manuf_cbp")
#' weight <- "timepwt48"
#' Z <- setdiff(names(ADH_local_wide), index)
#' G <- "trade_"
#'
#' bw(ADH_master, y, x, controls, weight,
#'    ADH_local_wide, Z, ADH_global, G,
#'    L_id = index, T_id = NULL)
#' @export
bw <- function(master, y, x, controls = NULL, weight = NULL,
              local, z, global, g, L_id, T_id = NULL) {
    # Parsing the master file
    Y <- as.matrix(master[y])
    X <- as.matrix(master[x])
    n <- nrow(Y)

    if (is.null(weight)) {
        weight <- as.matrix(rep(1, n))
    } else {
        weight <- as.matrix(master[weight])
    }

    if (is.null(controls)) {
        WW <- matrix(1, n, 1)
    } else {
        W <- as.matrix(master[controls])
        WW <- cbind(W, matrix(1, n, 1))
    }

    # Match the data, local and global matrices
    matched_matrices <- match_matrices(master, local, global, z, g, L_id, T_id)

    Z_matched <- matched_matrices[["Z"]]
    B <- matched_matrices[["B"]]
    Bk <- matched_matrices[["Bk"]]

    # Compute the coefficients by groups
    coeffs <- ComputeAlphaBeta(Y, X, WW, weight, Z_matched, B, Bk)

    # Return a tibble
    tibble::as_tibble(cbind(
        global,
        alpha = coeffs[["alpha"]], beta = coeffs[['beta']],
        gamma = coeffs[["gamma"]], pi = coeffs[['pi']])
    )
}

match_matrices <- function(data, local, global, z, g, L_id, T_id = NULL) {
    stopifnot(length(g) == 1 || !is.null(T_id))

    # Create Z and G matrices
    Z <- as.matrix(local[z])
    G <- as.matrix(global[g])

    # Create locations ids by observation
    data_L_index <- Reduce(paste0, data[L_id])
    rownames(Z) <- Reduce(paste0, local[L_id])

    # Identify the observations in the data by ids and periods
    data_T_index <- Reduce(paste0, data[T_id])
    if(ncol(G) == 1) data_LT_index <- paste0(data_L_index, colnames(G))
    else data_LT_index <- paste0(data_L_index, data_T_index)

    # In case these are numeric
    data_L_index <- as.character(data_L_index)
    data_LT_index <- as.character(data_LT_index)

    # Create all the possible ids L*T
    names_lt <- expand.grid(rownames(Z), colnames(G))
    names_lt <- paste0(names_lt[[1]], names_lt[[2]])

    # Match Z with data observations locations
    Z_matched <- Z[data_L_index, ]
    Z_matched <- as.matrix(Z_matched)

    # Row wise multiplication Z % G
    Bk_vec <- list()
    for (i in 1:ncol(G)) Bk_vec[[i]] <- G[, i] * t(Z)
    Bk_vec <- Reduce(cbind, Bk_vec) # dimension (K x L*T)
    colnames(Bk_vec) <- names_lt
    Bk_vec <- Bk_vec[, data_LT_index] # Reorder columns to match the data
    Bk_vec <- as.matrix(Bk_vec)

    # Matrix multiplication Z * G
    B_vec <- Z %*% G
    B_vec <- as.matrix(c(B_vec)) # dimension (L*T x 1)
    rownames(B_vec) <- names_lt
    B_vec <- B_vec[data_LT_index, ] # Reorder columns to match the data
    B_vec <- as.matrix(B_vec)

    list(B = B_vec, Bk = Bk_vec, Z = Z_matched)
}

ComputeAlphaBeta <- function(y, x, WW, weight, Z, B, Bk) {
    weightSQR <- sqrt(weight)

    for (i in 1:ncol(x)) x[, i] <- x[, i] * weightSQR
    for (i in 1:ncol(y)) y[, i] <- y[, i] * weightSQR
    for (i in 1:ncol(Z)) Z[, i] <- Z[, i] * weightSQR
    for (i in 1:ncol(WW)) WW[, i] <- WW[, i] * weightSQR
    for (i in 1:ncol(B)) B[, i] <- B[, i] * weightSQR
    for (i in 1:nrow(Bk)) Bk[i, ] <- Bk[i, ] * weightSQR
    rm(weightSQR)

    xx <- x - WW %*% qr.solve(WW, x)
    rm(x)
    yy <- y - WW %*% qr.solve(WW, y)
    rm(y)
    ZZ <- Z - WW %*% qr.solve(WW, Z)
    rm(WW)

    Zxx <- crossprod(Z, xx)
    Zyy <- crossprod(Z, yy)
    ZZZZ <- rowSums(t(ZZ) * t(ZZ))

    Alpha <- (Bk %*% xx) / as.numeric(crossprod(B, xx))
    Beta <- Zyy / Zxx
    Gamma <- Zyy / ZZZZ
    pi <- crossprod(ZZ, xx) / ZZZZ

    attributes(Alpha)[["dimnames"]] <- NULL
    attributes(Beta)[["dimnames"]] <- NULL
    attributes(Gamma)[["dimnames"]] <- NULL
    attributes(pi)[["dimnames"]] <- NULL

    list(alpha = Alpha, beta = Beta, gamma = Gamma, pi = pi)
}

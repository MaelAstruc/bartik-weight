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
#' @importFrom Rcpp sourceCpp
#' @useDynLib bartik.weight, .registration = TRUE
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
        alpha = coeffs[[1]], beta = coeffs[[2]],
        gamma = coeffs[[3]], pi = coeffs[[4]],
        se    = coeffs[[5]]
    ))
}

paste_ <- function(...) {
    paste(..., sep = "_")
}

match_matrices <- function(master, local, global, z, g, L_id, T_id = NULL) {
    stopifnot(length(g) == 1 || !is.null(T_id))

    # Create list of ids for the data
    master_l <- Reduce(paste_, master[L_id])
    master_t <- Reduce(paste_, master[T_id])
    master_lt <- Reduce(paste_, master[c(L_id, T_id)])

    # Create Z and G matrices
    Z <- as.matrix(local[z])
    G <- as.matrix(global[g])

    # Use row name to identify locations
    rownames(Z) <- Reduce(paste_, local[L_id])

    # Reduce Z and G
    Z <- Z[rownames(Z) %in% master_l, ]
    if (!is.null(T_id)) G <- G[, colnames(G) %in% master_t]

    # Create lists of id for the instrument
    names_l <- rownames(Z)
    names_k <- colnames(Z)

    if (is.null(T_id)) {
        names_lt <- names_l
        names_kt <- names_k
    } else {
        names_t <- colnames(G)
        names_lt <- expand.grid(names_l, names_t)
        names_lt <- Reduce(paste_, names_lt)
        names_kt <- expand.grid(names_k, names_t)
        names_kt <- Reduce(paste_, names_kt)
    }

    # In case of panel data Z is (L*T x T*K) with 0 when T does not match
    Z_T <- matrix(
        data = 0,
        nrow = nrow(Z) * ncol(G),
        ncol = ncol(Z) * ncol(G),
        dimnames = list(names_lt, names_kt)
    )

    for (t in 1:ncol(G)) {
        start_row <- (t - 1) * nrow(Z) + 1
        end_row <- t * nrow(Z)
        start_col <- (t - 1) * ncol(Z) + 1
        end_col <- t * ncol(Z)

        Z_T[start_row:end_row, start_col:end_col] <- Z
    }

    # Same for G
    G_T <- c(G)
    names(G_T) <- names_kt

    # Match Z with the data
    Z_matched <- Z_T[master_lt, ]

    # Normal Bartik
    B <- Z %*% G
    B <- as.matrix(c(B)) # dimension (L*T x 1)
    rownames(B) <- names_lt
    B <- B[master_lt, ] # Reorder rows to match the data
    B <- as.matrix(B)

    # Row wise multiplication Z % G
    Bk <- matrix(nrow = nrow(Z_T), ncol = ncol(Z_T), dimnames = dimnames(Z_T))

    for (kt in names(G_T)) Bk[, kt] <- G_T[kt] * Z_T[, kt]

    Bk <- t(Bk)
    Bk <- Bk[, master_lt] # Reorder columns to match the data

    return(list(B = B, Bk = Bk, Z = Z_matched))
}

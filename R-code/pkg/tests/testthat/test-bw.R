# Check that the matching in done right

test_that("match_matrices works", {
    # Simple data
    # N = 3, L = 2, K = 2, T = 2
    master <- data.frame(
        id = c(1, 1, 2, 2, 3, 3),
        location = paste0("location_", c(1, 1, 2, 1, 2, 2)),
        period = paste0("period_", c(1, 2, 1, 2, 1, 2))
    )

    local <- data.frame(
        location = paste0("location_", c(1, 2)),
        group_1 = c(1, 2),
        group_2 = c(3, 4)
    )

    global <- data.frame(
        group = paste0("group_", c(1, 2)),
        period_1 = c(1, 1),
        period_2 = c(2, 2)
    )

    z <- c("group_1", "group_2")
    g <- c("period_1", "period_2")

    L_id <- "location"
    T_id <- "period"

    result <- match_matrices(master, local, global, z, g, L_id, T_id)

    # Matrices multiplication results
    # B <- as.matrix(c(4, 6, 8, 12))
    # rownames(B) <- c("location_1period_1", "location_2period_1",
    #                  "location_1period_2", "location_2period_2")
    #
    # Bk <- as.matrix(data.frame(
    #     location_1period_1 = c(1, 3),
    #     location_2period_1 = c(2, 4),
    #     location_1period_2 = c(2, 6),
    #     location_2period_2 = c(4, 8)
    # ))
    # rownames(Bk) <- c("group_1", "group_2")

    # Matrices after matching
    Z_expected <- as.matrix(data.frame(
        group_1 = c(1, 1, 2, 1, 2, 2),
        group_2 = c(3, 3, 4, 3, 4, 4)
    ))
    rownames(Z_expected) <- paste0("location_", c(1, 1, 2, 1, 2, 2))

    B_expected <- as.matrix(c(4, 8, 6, 8, 6, 12))
    rownames(B_expected) <- c("location_1period_1", "location_1period_2",
                              "location_2period_1", "location_1period_2",
                              "location_2period_1", "location_2period_2")

    Bk_expected <- as.matrix(data.frame(
        location_1period_1 = c(1, 3),
        location_1period_2 = c(2, 6),
        location_2period_1 = c(2, 4),
        location_1period_2 = c(2, 6),
        location_2period_1 = c(2, 4),
        location_2period_2 = c(4, 8),
        check.names = FALSE
    ))
    rownames(Bk_expected) <- c("group_1", "group_2")

    expect_equal(result[["Z"]], Z_expected)
    expect_equal(result[["B"]], B_expected)
    expect_equal(result[["Bk"]], Bk_expected)
})

# Check that the coefficients are right
test_that("ComputeAlphaBeta works", {
    # Short data
    # N = 3, L = 2, K = 2, T = 2
    y <- as.matrix(c(1, 3, 2, 6, 3, 9))
    x <- as.matrix(c(2, 1, 3, 7, 8, 10))
    WW <- as.matrix(data.frame(
        fixed = rep(1, 6),
        control_1 = c(1, 2, 7, 2, 6, 4)
    ))
    weight <- as.matrix(1:6)

    # The results from the previous test
    Z <- as.matrix(data.frame(
        group_1 = c(1, 1, 2, 1, 2, 2),
        group_2 = c(3, 3, 4, 3, 4, 4)
    ))

    B <- as.matrix(c(4, 8, 6, 8, 6, 12))
    rownames(B) <- c("location_1period_1", "location_1period_2",
                              "location_2period_1", "location_1period_2",
                              "location_2period_1", "location_2period_2")

    Bk <- as.matrix(data.frame(
        location_1period_1 = c(1, 3),
        location_1period_2 = c(2, 6),
        location_2period_1 = c(2, 4),
        location_1period_2 = c(2, 6),
        location_2period_1 = c(2, 4),
        location_2period_2 = c(4, 8),
        check.names = FALSE
    ))

    result <- ComputeAlphaBeta(y, x, WW, weight, Z, B, Bk)

    # Expected results
    alpha <- as.matrix(c(0.39297336, 0.60702664))
    beta <- as.matrix(c(0.9588244, 0.9588244))
    gamma <- as.matrix(c(8.9140083 , 8.9140083))
    pi <- as.matrix(c(9.29681, 9.29681))

    expect_equal(result[[1]], alpha)
    expect_equal(result[[2]], beta)
    expect_equal(result[[3]], gamma)
    expect_equal(result[[4]], pi)
})

# Check that the alphas and betas are the same since @jjchern version

test_that("bw still gives the same output with ADH example", {
    bw_ADH_ref <- readRDS(testthat::test_path("test_data", "bw_ADH.rds"))
    ref_cols <- colnames(bw_ADH_ref)

    ADH_local %>%
        dplyr::mutate(ind = stringr::str_glue("t{year}_sh_ind_{ind}")) %>%
        tidyr::spread(ind, sh_ind_, fill = 0) -> ADH_local2

    # Prepare variables in the master tibble
    index = c("czone", "year")
    y = "d_sh_empl_mfg"
    x = "d_tradeusch_pw"

    controls = c("reg_midatl", "reg_encen", "reg_wncen", "reg_satl",
                 "reg_escen", "reg_wscen", "reg_mount", "reg_pacif", "l_sh_popedu_c",
                 "l_sh_popfborn", "l_sh_empl_f", "l_sh_routine33", "l_task_outsource",
                 "t2", "l_shind_manuf_cbp")

    weight = "timepwt48"

    # Prepare variables in the local tibble
    Z = setdiff(names(ADH_local_wide), index)

    # Prepare variables in the global tibble
    G = "trade_"

    # Estimate the weight (alpha) and the IV estimates (beta)
    bw_ADH_test = bw(ADH_master, y, x, controls, weight,
                     ADH_local2, Z, ADH_global, G,
                     L_id = index, T_id = NULL)

    expect_equal(bw_ADH_test[, ref_cols], bw_ADH_ref)
})


test_that("bw still gives the same output with BAR example", {
    bw_BAR_ref <- readRDS(testthat::test_path("test_data", "bw_BAR.rds"))
    ref_cols <- colnames(bw_BAR_ref)

    BAR_local %>%
        dplyr::select(-sh_ind_) %>%
        dplyr::mutate(ind = stringr::str_glue("t{year}_init_sh_ind_{ind}")) %>%
        tidyr::spread(ind, init_sh_ind_, fill = 0) -> BAR_local2

    # Prepare variables in the master tibble
    index = c("czone", "year")
    y = "wage_ch"
    x = "emp_ch"
    weight = "pop1980"
    controls = setdiff(names(BAR_master), c(index, y, x, weight))

    # Prepare variables in the local tibble
    Z = setdiff(names(BAR_local2), c(index))

    # Prepare variables in the global tibble
    G = "nat_empl_ind_"

    # Estimate the weight (alpha) and the IV estimates (beta)
    bw_BAR_test = bw(BAR_master, y, x, controls, weight,
                     BAR_local2, Z, BAR_global, G,
                     L_id = index, T_id = NULL)

    expect_equal(bw_BAR_test[, ref_cols], bw_BAR_ref)
})

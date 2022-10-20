# Check that the alphas and betas are the same since @jjchern version

test_that("bw still gives the same output with ADH example", {
    bw_ADH_ref <- readRDS(testthat::test_path("test_data", "bw_ADH.rds"))
    ref_cols <- colnames(bw_ADH_ref)

    ADH_local %>%
        dplyr::mutate(ind = stringr::str_glue("t{year}_sh_ind_{ind}")) %>%
        tidyr::spread(ind, sh_ind_, fill = 0) -> ADH_local2

    y = "d_sh_empl_mfg"
    x = "d_tradeusch_pw"

    controls = c("reg_midatl", "reg_encen", "reg_wncen", "reg_satl",
                 "reg_escen", "reg_wscen", "reg_mount", "reg_pacif", "l_sh_popedu_c",
                 "l_sh_popfborn", "l_sh_empl_f", "l_sh_routine33", "l_task_outsource",
                 "t2", "l_shind_manuf_cbp")

    weight = "timepwt48"

    # Prepare variables in the local tibble
    Z = setdiff(names(ADH_local_wide), c("czone", "year"))

    # Prepare variables in the global tibble
    G = "trade_"

    # Estimate the weight (alpha) and the IV estimates (beta)
    bw_ADH_test = bw(ADH_master, y, x, controls, weight, ADH_local2, Z, ADH_global, G)

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
    bw_BAR_test = bw(BAR_master, y, x, controls, weight, BAR_local2, Z, BAR_global, G)

    expect_equal(bw_BAR_test[, ref_cols], bw_BAR_ref)
})

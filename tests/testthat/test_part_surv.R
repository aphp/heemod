library(flexsurv)

surv_dist1 <- define_surv_dist(
  distribution = "exp",
  rate = .003
)
surv_dist2 <- define_surv_dist(
  distribution = "exp",
  rate = .002
)
surv_dist3 <- define_surv_dist(
  distribution = "exp",
  rate = .005
)
surv_dist4 <- define_surv_dist(
  distribution = "exp",
  rate = .004
)


suppressMessages(
  ps <- define_part_surv(
    pfs = join(surv_dist3, surv_dist4, at=365),
    os = join(surv_dist1, surv_dist2, at=365),
    cycle_length = c(365, 365)
  )
)
suppressMessages(
  ps1 <- define_part_surv(
    pfs = join(surv_dist3, surv_dist4, at=365) |> apply_hr(0.8),
    os = join(surv_dist1, surv_dist2, at=365) |> apply_hr(0.8),
    cycle_length = c(365, 365)
  )
)

sA <-  define_state(
  cost = 10000, ut = 1
)
sB <-  define_state(
  cost = 20000, ut = 1
)
sC <-  define_state(
  cost = 0, ut = 0
)

stratPS <- define_strategy(
  transition = ps,
  A = sA, B = sB, C = sC
)

stratPS1 <- define_strategy(
  transition = ps1,
  A = sA, B = sB, C = sC
)

suppressMessages(
  resPS <- run_model(
  Strat1 = stratPS,
  Strat2 = stratPS1,
  cycles = 10,
  cost = cost,
  effect = ut,
  method = "end"
))

test_that(
  "part surv works", {
    
    surv_dist_1 <- define_surv_dist(
      distribution = "exp", rate = 0.5)
    fit_cov <- flexsurvreg(
      formula = survival::Surv(rectime, censrec) ~ group,
      dist = "weibull", 
      data = bc) |> 
      define_surv_fit()
    fitcov_medium <- set_covariates(fit_cov, group = "Medium")
    km_cov <- survival::survfit(
      formula = survival::Surv(rectime, censrec) ~ group,
      data = bc) |> 
      define_surv_fit()
    km_medium <- set_covariates(km_cov, group = "Medium")
    
    suppressMessages({
      ps <- define_part_surv(
        pfs = surv_dist_1,
        os = km_medium |>
          join(fitcov_medium, 
                  at = 730),
        cycle_length = c(1, 365)  # 1 for pfs, 365 for os
      )})
    
    sA <- define_state(cost = 10, ut = 1)
    sB <- define_state(cost = 20, ut = 0.5)
    sC <- define_state(cost = 0, ut = 0)
    
    stratPS <- define_strategy(transition = ps, A = sA, B = sB, C = sC)
    
    param <- define_parameters(
      p1 = compute_surv(
        surv_dist_1,
        time = model_time  # can also be state_time
      ),
      p2 = compute_surv(
        km_medium |>
          join(fitcov_medium, 
                  at = 730),
        time = model_time, cycle_length = 365  # time is in days in km_medium, in years in model_time
      ))
    
    tm <- define_transition(
      1 - p1, C, p2,
      0, C, p2,
      0, 0, 1)
    #> No named state -> generating names.
    
    stratTM <- define_strategy(
      transition = tm, A = sA, B = sB, C = sC)
    
    suppressWarnings({
      resPS <- run_model(
        stratPS,
        cycles = 10)
    })
    suppressWarnings({
      resTM <- run_model(
        parameters = param,
        stratTM,
        cycles = 10)
    })
    
    expect_equal(
      get_counts(resPS)$count,
      get_counts(resTM)$count
    )
    
    suppressWarnings({
      resPS <- run_model(
        stratPS,
        cycles = 10,
        method = "end")
    })
    suppressWarnings({
      resTM <- run_model(
        parameters = param,
        stratTM,
        cycles = 10,
        method = "end")
    })
    
    expect_equal(
      get_counts(resPS)$count,
      get_counts(resTM)$count
    )
  }
)

suppressMessages({
  suppressMessages(
  ps <- define_part_surv(
    pfs = join(surv_dist3),
    os = join(surv_dist1),
    state_names = c("ProgressionFree", "Progressive", "Death"),
    cycle_length = c(365, 365)
  ))
})
test_that(
  "errors with inappropriate state names", {
    expect_error(
      suppressMessages(
      define_part_surv(
        pfs = join(surv_dist3),
        os = join(surv_dist1),
        state_names = c("NoDisease", "Progressive", "Death"),
        cycle_length = c(365, 365)
      )),
      "Progression free state (only) must have 'free' in its name",
      fixed = TRUE
    )
    expect_error(
      suppressMessages(
      define_part_surv(
        pfs = join(surv_dist3),
        os = join(surv_dist1),
        state_names = c("ProgressionFree", "Progressive", "Kaput"),
        cycle_length = c(365, 365)
      )),
      "State name representing death"
    )
    expect_error(
      suppressMessages(
      define_part_surv(
        pfs = join(surv_dist3),
        os = join(surv_dist1),
        state_names = c("ProgressionFree", "Progressive",
                        "uh-oh", "Death"),
        cycle_length = c(365, 365)
      )),
      "If there are 4 states, a state must be called 'terminal'",
      fixed = TRUE
    )
    expect_error(
      suppressMessages(
      define_part_surv(
        pfs = join(surv_dist3),
        os = join(surv_dist1),
        state_names = c(
          "ProgressionFree",
          "Progressivebutfree",
          "terminal",
          "Death"
        ),
        cycle_length = c(365, 365)
      )),
      "Progression free state (only) must have 'free' in its name",
      fixed = TRUE
    )
  })

test_that(
  "join_fits_to_def works", {
    
    surv_def <- read_file(system.file("tabular/surv", 
                                      "use_fits.csv", 
                                      package = "heemod"))
    surv_def$.subset <- "all"
    fake_fit_tib <- read_file(system.file("tabular/surv",
                                          "fake_fit_tib.csv", 
                                          package = "heemod"))
    state_names <- c("ProgressionFree", "ProgressiveDisease", 
                     "Terminal", "Death")
    ## basically just make sure it runs, since we're using fake fits
    zz <- join_fits_to_def(surv_def, fake_fit_tib)
    expect_identical(names(zz), c(".strategy", ".type", "dist",
                                  ".subset", "fit", "set_def"))
    surv_def_join <- read_file(system.file("tabular/surv", 
                                           "use_fits_join.csv", 
                                           package = "heemod"))
    bad_surv_def <- surv_def_join
    bad_surv_def[[1, "dist"]] <- "fit('bad')"
    expect_error(capture.output(join_fits_to_def(
      bad_surv_def, 
      fake_fit_tib)),
      "disallowed distribution name in use_fits")
    expect_error(capture.output(join_fits_to_def(
      surv_def_join,
      fake_fit_tib)),
      "subset name GT50 from specification of fits"
    )
    names(surv_def)[1] <- "strategy"
    expect_error(join_fits_to_def(surv_def, fake_fit_tib),
                 "missing required names in 'surv_def':", fixed = TRUE)
    names(surv_def)[1] <- ".strategy"
    names(fake_fit_tib)[1] <- ".strategy"
    expect_error(join_fits_to_def(surv_def, fake_fit_tib),
                 "missing required names in 'fit_tibble':", fixed = TRUE)
    names(fake_fit_tib)[1] <- "type"
    fake_fit_tib$set_name[c(15)] <- "GT50"
    expect_error(capture.output(join_fits_to_def(surv_def_join, 
                                  fake_fit_tib)),
                 "fit not found for line"
    )
      }
  )

test_that("we catch bad names in construct_part_surv_tib",
          {
            surv_def <- read_file(system.file("tabular/surv", 
                                              "use_fits.csv", 
                                              package = "heemod"))
            surv_def$.subset <- "all"
            state_names <- c("ProgressionFree", "ProgressiveDisease", 
                             "Terminal", "Death")

            names(surv_def)[1] <- "strategy"
            expect_error(construct_part_surv_tib(surv_def, NULL, state_names),
                         "missing required names in 'surv_def':", fixed = TRUE)
            names(surv_def)[1] <- ".strategy"
            names(surv_def)[2] <- "type"
            expect_error(construct_part_surv_tib(surv_def, NULL, state_names),
                         "missing required names in 'surv_def':", fixed = TRUE)
            names(surv_def)[2] <- ".type"
            names(surv_def)[3] <- "DIST"
            expect_error(construct_part_surv_tib(surv_def, NULL, state_names),
                         "missing required names in 'surv_def':", fixed = TRUE)
            names(surv_def)[3] <- "dist"
          }
)

test_that("we can run construct_part_surv_tib",
          {
            use_fits <- read_file(system.file("tabular/surv",
                                  "example_use_fits_explicit_dists.csv",
                                  package = "heemod"))
            ref <- read_file(system.file("tabular/surv",
                                         "example_oncSpecs_explicit_dists.csv",
                                         package = "heemod"))
            suppressMessages(explicit_dist_part_surv <- 
              construct_part_surv_tib(use_fits, ref,             
                                    state_names <- c("ProgressionFree", 
                                                     "ProgressiveDisease", 
                                                      "Terminal", "Death")
            ))
            for_A <- dplyr::filter(explicit_dist_part_surv, .strategy == "A")
            for_B <- dplyr::filter(explicit_dist_part_surv, .strategy == "B")
            expect_equal(round(compute_surv(for_A[["part_surv"]][[1]]$pfs, 1), 4), 0.0332)
            expect_equal(round(compute_surv(for_A[["part_surv"]][[1]]$os, 1), 4), 0.0150)
            expect_equal(round(compute_surv(for_B[["part_surv"]][[1]]$pfs, 1), 4), 0.0472)
            expect_equal(round(compute_surv(for_B[["part_surv"]][[1]]$os, 1), 4), 0.0213)

            use_fits <- read_file(system.file("tabular/surv",
                                              "use_fits_mixed.csv",
                                              package = "heemod"))
            ref <- read_file(system.file("tabular/surv",
                                         "example_oncSpecs_mixed.csv",
                                         package = "heemod"))
            ref$full_file <- file.path(system.file("tabular/surv", package = "heemod"),
                                                   ref$file)
            suppressMessages( mixed_dist_part_surv <- 
              construct_part_surv_tib(use_fits, ref,             
                                      state_names <- c("ProgressionFree", 
                                                       "ProgressiveDisease", 
                                                       "Terminal", "Death")
              ))
            expect_identical(class(mixed_dist_part_surv[["part_surv"]][[1]]$os)[1],
                             "quosure")
            expect_identical(eval_tidy(mixed_dist_part_surv[["part_surv"]][[1]]$pfs),
                             'define_surv_dist(distribution = "exp", rate = 1/100)')
            expect_identical(class(eval_tidy(mixed_dist_part_surv[["part_surv"]][[1]]$os)),
                             "flexsurvreg")
            prob <- compute_surv(eval_tidy(mixed_dist_part_surv[["part_surv"]][[1]]$os), 1)
            expect_equal(round(prob, 5), 0.00213)
                      })

test_that("join_fits_across_time works", 
          {
    surv_def_join <- read_file(system.file("tabular/surv", 
                                           "use_fits_join.csv", 
                                           package = "heemod"))
    surv_def_join$fit <- letters[1:nrow(surv_def_join)]
    zz <- join_fits_across_time(surv_def_join[1:2,])
    surv_def_join <- dplyr::filter(surv_def_join, .subset == "all")
    zz <- join_fits_across_time(surv_def_join[1:2,])
    expect_error(capture.output(join_fits_across_time(
      surv_def_join[1:2,1:3])),
      "unless 'until' is also specified", fixed = TRUE)
  }
)


## don't have survival_from_data here.  Need to get saved data to run this test.

test_that("making part_surv from survival fits works",
          {
            location <- system.file("tabular/surv", package = "heemod")
            ok_surv_info <-
              heemod:::read_file(system.file("tabular/surv/survival_info.csv",
                                             package = "heemod"))
            fake_fit_tib <- read_file(system.file("tabular/surv",
                                                  "fake_fit_tib.csv", 
                                                  package = "heemod"))
            
            suppressMessages( these_part_survs <-
              part_survs_from_surv_inputs(fake_fit_tib,
                                          c("ProgressionFree",
                                            "Progressive",
                                            "Terminal",
                                            "Death")))
            expect_equal(nrow(these_part_survs), 14)
            expect_identical(
              names(these_part_survs),
              c("treatment", "set_name", "dist", "set_def", "part_surv")
            )
            expect_identical(class(these_part_survs$part_surv[[1]]), "part_surv")
          }
          )

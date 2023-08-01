library(survival)

surv_dist_1 <- define_surv_dist(
  distribution = "exp",
  rate = .4
)
surv_dist_2 <- define_surv_dist(
  distribution = "weibullPH",
  shape = 1
)

km_1 <- survfit(
  formula = Surv(time, status) ~ rx,
  data = colon) |>
  define_surv_fit()


test_that("define_psa gives correct output", {
  psa <- define_psa(surv_dist_1 ~ resample_surv(1000))
  expect_equal(psa$list_qdist$surv_dist_1, quote(resample_surv(1000)), ignore_attr = TRUE)
  expect_equal(psa$surv, "surv_dist_1", ignore_attr = TRUE)
  psa <- define_psa(surv_dist_1 ~  resample_surv(1000),
                    g + h ~ multinomial(5,12))
  expect_equal(psa$list_qdist$surv_dist_1, quote(resample_surv(1000)), ignore_attr = TRUE)
  expect_equal(names(psa$list_qdist), c("surv_dist_1", "g", "h"))
  expect_equal(psa$surv, "surv_dist_1", ignore_attr = TRUE)
})

test_that("psa surv_dist is correct", {

  param <- define_parameters(
    p1 = heemod:::compute_surv_(
      surv_dist_1,
      time = model_time
    ),
    p2 =  heemod:::compute_surv_(
      surv_dist_2,
      time = model_time
    )
  )

  tm <- define_transition(
    C, p1,
    0, C
  )

  tm2 <- define_transition(
    C, p2,
    0, C
  )

  sA <-  define_state(
    cost = 10, ut = 1
  )
  sB <-  define_state(
    cost = 20, ut = .5
  )

  stratTM <- define_strategy(
    transition = tm,
    A = sA, B = sB
  )

  stratTM2 <- define_strategy(
    transition = tm2,
    A = sA, B = sB
  )

  resTM <- run_model(
    parameters = param,
    stratTM,
    cycles = 15,
    cost = cost, effect = ut
  )

  psa <- define_psa(surv_dist_1 ~
                      resample_surv(1000))
  psaweak <- define_psa(surv_dist_1 ~
                          resample_surv(10))

  set.seed(1)
  resPSA <- run_psa(resTM, psa, 10)
  expect_true(
    abs(resTM$run_model$.cost - resPSA$run_model$cost) < 500
  )
  resPSAweak <- run_psa(resTM, psaweak, 10)
  expect_true(
    sd(resPSAweak$psa$cost)/sd(resPSA$psa$cost) > 5
  )
#   
#   # test_that("join_surv works with run_psa with use_surv_dist and resample_surv", {
#   #   join_surv <- join(surv_dist_1, surv_dist_2, surv_dist_1, at = c(2,5))
#   #   
#   #   param <- define_parameters(
#   #     p1 = heemod:::compute_surv_(
#   #       join_surv,
#   #       time = model_time
#   #     )
#   #   )
#   #   
#   #   resTM <- run_model(
#   #     parameters = param,
#   #     stratTM,
#   #     cycles = 15,
#   #     cost = cost, effect = ut
#   #   )
#   #   set.seed(1)
#   #   psa <- define_psa(surv_dist_1 ~ use_surv_dist(rexp(10000, 0.1)),
#   #                     surv_dist_2 ~ use_surv_dist(rexp(10000, 1)))
#   #   
#   #   
#   #   resPSA <- run_psa(resTM, psa, 10, TRUE)
#   #   for (i in 1:10){
#   #     pars <- resPSA$full[[1]][[i]]$parameters$p1
#   #     expect_true(all(pars[c(1,2, 6:15)] < 0.5))
#   #     expect_true(all(pars[c(3:5)] > 0.5))
#   #   }
#   # })
#   # set.seed(1)
#   # psa <- define_psa(surv_dist_1 ~ resample_surv(10000),
#   #                   surv_dist_2 ~ resample_surv(10000))
#   # 
#   # 
#   # resPSA <- run_psa(resTM, psa, 10, TRUE)
#   # for (i in 1:10){
#   #   pars <- resPSA$full[[1]][[i]]$parameters$p1
#   #   expect_true(all(pars[c(1,2, 6:15)] < 0.5))
#   #   expect_true(all(pars[c(3:5)] > 0.4))
#   # }
#   # 
#   # 
#   # test_that("dots is working for join_surv",{
#   #   
#   #   join_surv2 <- join_surv |>
#   #     join(surv_dist_2, at = 10)
#   #   
#   #   param <- define_parameters(
#   #     p1 = heemod:::compute_surv_(
#   #       join_surv2,
#   #       time = model_time
#   #     )
#   #   )
#   #   
#   #   resTM <- run_model(
#   #     parameters = param,
#   #     stratTM,
#   #     cycles = 15,
#   #     cost = cost, effect = ut
#   #   )
#   #   
#   #   resPSA <- run_psa(resTM, psa, 10, TRUE)
#   #   pars <- resPSA$full[[1]][[1]]$parameters$p1
#   #   expect_true(all(pars[c(1,2, 6:9)] < 0.5))
#   #   expect_true(all(pars[c(3:5)] > 0.4))
#   #   expect_true(all(pars[c(11:15)] == 0))
#   # })
})

test_that("psa surv_fit is correct", {
  param <- define_parameters(
    p1 = heemod:::compute_surv_(
      km_1,
      time = model_time,
      cycle_length = 30
    )
  )

  tm <- define_transition(
    C, p1,
    0, C
  )

  tm2 <- define_transition(
    C, p2,
    0, C
  )

  sA <-  define_state(
    cost = 10, ut = 1
  )
  sB <-  define_state(
    cost = 20, ut = .5
  )

  stratTM <- define_strategy(
    transition = tm,
    A = sA, B = sB
  )
  resTM <- run_model(
    parameters = param,
    stratTM,
    cycles = 15,
    cost = cost, effect = ut
  )
  set.seed(1)
  psa <- define_psa(km_1 ~ resample_surv())
  resPSA <- run_psa(resTM, psa, 10)
  expect_length(unique(resPSA$psa$ut), 10)

  expect_gt(sd(resPSA$psa$cost), 0)
  expect_gt(sd(resPSA$psa$ut), 0)
  expect_equal(identical(resPSA$run_model, resPSA$model$run_model), FALSE)
  
  expect_true(
    sd(resPSA$psa$cost) > 0
  )

   km_2 <- km_1 |>
     set_covariates(rx = "Lev+5FU")
   
   psa <- define_psa(km_2 ~ resample_surv())
   
   expect_error(suppressWarnings(run_psa(resTM, psa, 10)), "No defined parameter")
   
   psa <- define_psa(km_2 ~ resample_surv(),
                     km_1 ~ resample_surv()
                     )
   
   expect_warning(run_psa(resTM, psa, 10), "km_2 not previously defined")
  ### do some tests here
})
# 


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
  expect_equal(psa$list_qdist$surv_dist_1, quote(resample_surv_dist(1000)), ignore_attr = TRUE)
  expect_equal(psa$surv, "surv_dist_1", ignore_attr = TRUE)
  psa <- define_psa(surv_dist_1 ~  resample_surv(1000),
                    g + h ~ multinomial(5,12))
  expect_equal(psa$list_qdist$surv_dist_1, quote(resample_surv_dist(1000)), ignore_attr = TRUE)
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
  
  test_that("survival operations work with run_psa and surv_dist", {
    modified_surv <- apply_or(surv_dist_1,3)
    param <- define_parameters(
          p1 = heemod:::compute_surv_(
            modified_surv,
            time = model_time
          ), 
          p2 = heemod:::compute_surv_(
            surv_dist_1,
            time = model_time
          )
        )

        resTM <- run_model(
          parameters = param,
          stratTM2,
          stratTM,
          cycles = 15,
          cost = cost, effect = ut
        )

        psa <- define_psa(modified_surv ~ resample_surv(10000))
      
        expect_error(run_psa(resTM, psa, 10), "the initial survival distribution")
    
        psa <- define_psa(surv_dist_1 ~ resample_surv(10000))
        
        set.seed(1)            
        resPSA <- run_psa(resTM, psa, 10)
        
        tbl <- resPSA$psa %>%
          group_by(.strategy_names) %>%
          dplyr::summarise(m = mean(.cost),
                    sd = sd(.cost))
        
        expect_true(all(tbl$sd > 0))
        
        expect_true(all(abs(resPSA$run_model$cost - resPSA$model$run_model$cost) < 100))
        
        
        set.seed(1)
        list_new <- eval_resample(psa, 10)$surv_dist_1
        list_surv_1 <- map(list_new, function(x){
          structure(x[[1]], class = c("surv_dist", "surv_object"))
        })
        list_surv_modified <- map(list_new, function(x){
          structure(list(dist = x[[1]], or = 3), class = c("surv_po", "surv_object"))
        })
       res <- purrr::map2(list_surv_1, list_surv_modified, function(x, y){
          surv_dist_1 <- x
          modified_surv <- y
          run_model(
            parameters = param,
            stratTM2,
            stratTM,
            cycles = 15,
            cost = cost, effect = ut
          ) %>% 
            `[[`("run_model")
        }) %>% 
         dplyr::bind_rows() %>% 
         dplyr::arrange(.strategy_names) %>% 
         dplyr::select(cost, ut) %>% 
         as_tibble()
       expect_equal(resPSA$psa %>% 
         dplyr::select(cost, ut), res)

  })

test_that("join_surv works with run_psa with use_surv_dist and resample_surv", {
  join_surv <- join(surv_dist_1, surv_dist_2, surv_dist_1, at = c(2,5))

  param <- define_parameters(
    p1 = heemod:::compute_surv_(
      join_surv,
      time = model_time
    )
  )

  resTM <- run_model(
    parameters = param,
    stratTM,
    cycles = 15,
    cost = cost, effect = ut
  )
  set.seed(1)
  psa <- define_psa(surv_dist_1 ~ resample_surv(10000),
                    surv_dist_2 ~ resample_surv(1000))


  resPSA <- run_psa(resTM, psa, 10, keep = TRUE)
  set.seed(1)
  for (i in 1:10){
    pars <- resPSA$full[[1]][[i]]$parameters$p1
    expect_true(all(pars[c(1,2, 6:15)] < 0.5))
    expect_true(all(pars[c(3:5)] > 0.4))
  }
})


  test_that("dots is working for join_surv",{

    join_surv2 <- join_surv |>
      join(surv_dist_2, at = 10)

    param <- define_parameters(
      p1 = heemod:::compute_surv_(
        join_surv2,
        time = model_time
      )
    )

    resTM <- run_model(
      parameters = param,
      stratTM,
      cycles = 15,
      cost = cost, effect = ut
    )

    resPSA <- run_psa(resTM, psa, 10, keep = TRUE)
    pars <- resPSA$full[[1]][[1]]$parameters$p1
    expect_true(all(pars[c(1,2, 6:9)] < 0.5))
    expect_true(all(pars[c(3:5)] > 0.4))
    expect_true(all(pars[c(11:15)] > 0.4))
  })
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

  
  test_that("errors run_psa resample", {
   km_2 <- km_1 |>
     set_covariates(rx = "Lev+5FU")
   
   psa <- define_psa(km_2 ~ resample_surv())
   
   expect_error(suppressWarnings(run_psa(resTM, psa, 10)), "No defined parameter")
   
   psa <- define_psa(km_2 ~ resample_surv(),
                     km_1 ~ resample_surv()
                     )
   
   expect_warning(run_psa(resTM, psa, 10), "km_2 neither used by define_parameters")
   
   param <- define_parameters(
     p1 = heemod:::compute_surv_(
       km_2,
       time = model_time,
       cycle_length = 30
     )
   )
   
   resTM <- run_model(
     parameters = param,
     stratTM,
     cycles = 15,
     cost = cost, effect = ut
   )
   expect_error(run_psa(resTM, psa, 10), 
                cli::cli_text("`km_2` must be a {.cls surv_fit}"))
   
   psa <- define_psa(km_1 ~ resample_surv(n = 100))
   expect_warning(eval_resample(psa, 2), "should not contain the `n` argument")
   
    survfit(
     formula = Surv(time, status) ~ 1,
     data = survival::colon) |>
     define_surv_fit() |>
    expect_message("include the package environment")
  })
  
  test_that("survival operations work with run_psa and surv_fit", {
    modified_surv <- apply_or(km_1,3) |>
      set_covariates(rx = "Lev+5FU")
    
    param <- define_parameters(
      p1 = heemod:::compute_surv_(
        km_1,
        time = model_time,
        cycle_length = 30
      ),
        p2 = heemod:::compute_surv_(
          modified_surv,
          time = model_time,
          cycle_length = 30
      )
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
      stratTM2,
      stratTM,
      cycles = 15,
      cost = cost, effect = ut
    )
    
    psa <- define_psa(km_1 ~ resample_surv())
    
    set.seed(1)            
    resPSA <- run_psa(resTM, psa, 10)
    
    tbl <- resPSA$psa %>%
      group_by(.strategy_names) %>%
      dplyr::summarise(m = mean(.cost),
                       sd = sd(.cost))
    
    expect_true(all(tbl$sd > 0))
    
    expect_true(all(abs(resPSA$run_model$cost - resPSA$model$run_model$cost) < 5000))
    
    
    set.seed(1)
    list_new <- eval_resample(psa, 10)$km_1
    colons <- map(list_new, function(x){
      eval_tidy(x[[1]])
    })
    res <- purrr::map(colons, function(x){
      assign("colon", x, getOption("heemod.env"))
      mod <- run_model(
        parameters = param,
        stratTM2,
        stratTM,
        cycles = 15,
        cost = cost, effect = ut
      ) %>% 
        `[[`("run_model")
      rm(colon, envir = getOption("heemod.env"))
      mod
    }) %>% 
      dplyr::bind_rows() %>% 
      dplyr::arrange(.strategy_names) %>% 
      dplyr::select(cost, ut) %>% 
      as_tibble()
    expect_equal(resPSA$psa %>% 
                   dplyr::select(cost, ut), res)
    
  })
  
  test_that("modifiers work with a variable name", {
    tm <- define_transition(
      C, p1,
      0, C
    )
    modified_surv <- apply_hr(km_1, hr) 
    
    param <- define_parameters(
      hr = 3,
      p1 = heemod:::compute_surv_(
        modified_surv,
        time = model_time,
        cycle_length = 30
      )
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
    psa <- define_psa(hr ~ lognormal(3, 1))
    
    set.seed(1)            
    resPSA <- run_psa(resTM, psa, 10)
    expect_gt(sd(resPSA$psa$cost), 0)
  })
})
# 


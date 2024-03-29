test_that(
  "Model definition", {
    mat1 <- define_transition(
      state_names = c("X1", "X2"),
      1-a, a,
      1-b, b
    )
    s1 <- define_state(
      x = 234,
      y = 123
    )
    s2 <- define_state(
      x = 987,
      y = 1726
    )
    mod1 <- define_strategy(
      transition = mat1,
      X1 = s1,
      X2 = s2
    )
    mat2 <- define_transition(
      1-a, a,
      1-b, b
    )
    mod2 <- define_strategy(
      transition = mat2,
      s1,
      s2
    )
    expect_output(
      print(mod1),
      "A Markov model strategy:

    2 states,
    2 state values",
      fixed = TRUE
    )
    expect_output(
      str(mod1),
      "List of 3
 $ transition",
      fixed = TRUE
    )
    expect_output(
      print(names(mod2$states)),
      '"A" "B"',
      fixed = TRUE
    )
    expect_error(
      define_strategy(
        transition = mat2,
        s1,
        s2,
        s2
      )
    )
    expect_error(
      define_strategy(
        transition = mat1,
        X1 = s1,
        X3 = s2
      )
    )
  }
)

test_that(
  "Model evaluation, 1 model", {
    par1 <- define_parameters(
      a = .1,
      b = 1 / (model_time + 1)
    )
    mat1 <- define_transition(
      state_names = c("X1", "X2"),
      1-a, a,
      1-b, b
    )
    s1 <- define_state(
      x = 234,
      y = 123
    )
    s2 <- define_state(
      x = 987,
      y = 1726
    )
    mod1 <- define_strategy(
      transition = mat1,
      X1 = s1,
      X2 = s2
    )
    e_mod <- run_model(
      mod1,
      parameters = par1,
      init = c(1, 0),
      cycles = 5,
      cost = x,
      effect = y,
      method = "end"
    )
    expect_equal(
      round(e_mod$run_model$x), 1593
    )
    expect_equal(
      round(e_mod$run_model$y), 1515
    )
    expect_equal(
      e_mod$run_model$.strategy_names, "I"
    )
    expect_equal(
      round(e_mod$run_model$.cost), 1593
    )
    expect_equal(
      round(e_mod$run_model$.effect), 1515
    )
    expect_equal(
      heemod:::get_eval_init(heemod:::get_eval_strategy_list(e_mod)[[1]]),
      c(1, 0),
      ignore_attr = TRUE
    )
    
    s_mod <- summary(e_mod)
    expect_equal(
      round(s_mod$res_values$x), 1593
    )
    expect_equal(
      round(s_mod$res_values$y), 1515
    )
    expect_equal(
      s_mod$res_values$.strategy_names, "I"
    )
    
    expect_equal(
      nrow(s_mod$res_comp), 1
    )
    
    expect_error(
      run_model(
        mod1,
        init = c(1, 0, 0),
        cycles = 5
      )
    )
    expect_error(
      run_model(
        mod1,
        init = c(-1, 0),
        cycles = 5
      )
    )
    expect_error(
      run_model(
        mod1,
        init = c(-1, 0),
        cycles = -5
      )
    )
  }
)


test_that(
  "Model evaluation, 2 models", {
    
    par1 <- define_parameters(
      a = .1,
      b = 1 / (model_time + 1)
    )
    mat1 <- define_transition(
      state_names = c("X1", "X2"),
      1-a, a,
      1-b, b
    )
    s1 <- define_state(
      x = 234,
      y = 123
    )
    s2 <- define_state(
      x = 987,
      y = 1726
    )
    mod1 <- define_strategy(
      transition = mat1,
      X1 = s1,
      X2 = s2
    )
    mod2 <- define_strategy(
      transition = mat1,
      X1 = s1,
      X2 = s1
    )
    
    e_mod2 <- run_model(
      mod1, mod2,
      parameters = par1,
      init = c(1, 0),
      cycles = 5,
      cost = x,
      effect = y,
      method = "end"
    )
    
    expect_equal(
      round(e_mod2$run_model$x), c(1593, 1170)
    )
    expect_equal(
      round(e_mod2$run_model$y), c(1515, 615)
    )
    expect_equal(
      e_mod2$run_model$.strategy_names, c("I", "II")
    )
    expect_equal(
      round(e_mod2$run_model$.cost), c(1593, 1170)
    )
    expect_equal(
      round(e_mod2$run_model$.effect), c(1515, 615)
    )
    expect_equal(
      heemod:::get_eval_init(heemod:::get_eval_strategy_list(e_mod2)[[1]]),
      c(1, 0),
      ignore_attr = TRUE
    )
    
    s_mod2 <- summary(e_mod2)
    expect_equal(
      round(s_mod2$res_values$x),  c(1593, 1170)
    )
    expect_equal(
      round(s_mod2$res_values$y), c(1515, 615)
    )
    expect_equal(
      s_mod2$res_values$.strategy_names, c("I", "II")
    )
    
    expect_equal(
      nrow(s_mod2$res_comp), 2
    )
  }
)

test_that(
  "eval_matrix works", {
    par <- tibble::tibble(
      model_time = 2:3,
      a = c(.1, .2)
    )
    mat <- define_transition(
      C, 1/model_time,
      a, 1-a
    )
    
    res <- heemod:::eval_transition(mat, par)
    
    expect_identical(
      round(res[[1]], 2),
      structure(c(0.5, 0.1, 0.5, 0.9), .Dim = c(2L, 2L))
    )
    expect_identical(
      round(res[[2]], 2),
      structure(c(0.67, 0.2, 0.33, 0.8), .Dim = c(2L, 2L))
    )
    
    mat2 <- define_transition(
      C, C,
      a, 1-a
    )
    expect_error(
      heemod:::eval_matrix(mat2, par)
    )
  }
)

test_that(
  "compute_counts fails when needed", {
    lm <- structure(list(
      structure(c(0.5, 0.1, 0.5, 0.9),
                .Dim = c(2L, 2L)),
      structure(c(0.67, 0.2, 0.33, 0.8),
                .Dim = c(2L, 2L))),
      class = c("eval_matrix", "list"),
      state_names = c("A", "B"))
    
    expect_error(
      heemod:::compute_counts(
        lm, init = c(10, 0, 0), method = "beginning")
    )
    expect_error(
      heemod:::compute_counts(
        lm, init = c(10), method = "beginning")
    )
    expect_error(
      heemod:::compute_counts(
        lm, init = c(10, 0), method = "endzzz")
    )
  }
)

test_that(
  "compute_counts works", {
    lm <- structure(list(
      structure(c(0.5, 0.1, 0.5, 0.9),
                .Dim = c(2L, 2L)),
      structure(c(0.67, 0.2, 0.33, 0.8),
                .Dim = c(2L, 2L))),
      class = c("eval_matrix", "list"),
      state_names = c("A", "B"))
    infw <- data.frame(c(0, 0), c(0, 0))
    
    expect_equal(names(heemod:::compute_counts(
      lm, init = c(10, 0), method = "end", inflow = infw)), c("counts", "diff"))
    expect_identical(
      dim(heemod:::compute_counts(
        lm, init = c(10, 0), method = "beginning", inflow = infw) |> 
          heemod:::correct_counts(method = "beginning") |>
          magrittr::extract2("counts")),
      c(2L, 2L)
    )
    expect_identical(
      dim(heemod:::compute_counts(
        lm, init = c(10, 0), method = "end", inflow = infw) |> 
          heemod:::correct_counts(method = "end") |>
          magrittr::extract2("counts")) ,
      c(2L, 2L)
    )
    expect_identical(
      dim(heemod:::compute_counts(
        lm, init = c(10, 0), method = "life-table", inflow = infw) |> 
          heemod:::correct_counts(method = "life-table") |>
          magrittr::extract2("counts")),
      c(2L, 2L),
      ignore_attr = TRUE
    )
    
    expect_equal(
      unlist(heemod:::compute_counts(
        lm, init = c(10, 0), method = "beginning", inflow = infw) |> 
          heemod:::correct_counts(method = "beginning") |>
          magrittr::extract2("counts")),
      c(10, 5, 0, 5),
      ignore_attr = TRUE
    )
    expect_equal(
      unlist(heemod:::compute_counts(
        lm, init = c(10, 0), method = "end", inflow = infw) |> 
          heemod:::correct_counts(method = "end") |>
          magrittr::extract2("counts")),
      c(5.00, 4.35, 5.00, 5.65),
      ignore_attr = TRUE
    )
    expect_equal(
      unlist(heemod:::compute_counts(
        lm, init = c(10, 0), method = "life-table", inflow = infw) |> 
          heemod:::correct_counts(method = "life-table") |>
          magrittr::extract2("counts")),
      c(7.500, 4.675, 2.500, 5.325),
      ignore_attr = TRUE
    )
  }
)

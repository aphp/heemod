test_that(
  "State definition", {
    s1 <- define_state(
      x = 234,
      y = 123
    )
    s2 <- define_state(
      x = 987,
      y = 1726
    )
    expect_equal(names(s1), c(".dots", "starting_values"))
    expect_equal(names(s1), names(s2))
    expect_output(
      str(s1$.dots),
      'List of 2
 $ x: language ~234',
      fixed = TRUE
    )
    expect_output(
      print(s1),
      'A state with 2 values.

x = 234
y = 123',
      fixed = TRUE
    )
    expect_error(
      define_state(234, 54)
    )
    expect_output(
      str(
        modify(
          s1,
          x = 111
        )
      ),
      "List of 2
 $ x: language ~111",
      fixed = TRUE
    )
    expect_error(
      modify(
        s1,
        z = 111
      )
    )
    expect_error(
      define_state(
        model_time = 876
      )
    )
    expect_error(
      modify(
        s1,
        model_time = 678
      )
    )
    expect_equal(
      get_state_value_names(s1),
      c("x", "y")
    )
    expect_output(print(s2 |>
                          modify(starting_values = define_starting_values(x = 100))),
                  "A state with 2 values and 1 starting value.

x = 987
y = 1726
Start
x = 100", fixed = TRUE)
    expect_error(s2 |>
                   modify(starting_values = list(x = 100)), "Incorrect length", fixed = TRUE)
  }
)


test_that(
  "State list definition", {
    s1 <- define_state(
      x = 234,
      y = 123
    )
    s2 <- define_state(
      x = 987,
      y = 1726
    )
    sl1 <- heemod:::define_state_list(
      X1 = s1,
      X2 = s2
    )
    sl2 <- heemod:::define_state_list(
      s1,
      s2
    )
    expect_output(
      str(sl1),
      "List of 2
 $ X1:List of 2
  ..$ .dots          :List of 2
  .. ..$ x: language ~234",
      fixed = TRUE
    )
    expect_output(
      print(sl1),
      "A list of 2 states with 2 values each.

State names:

X1
X2

State values:

x
y"
    )
    expect_output(
      str(sl2),
      "List of 2
 $ A:List of 2
  ..$ .dots          :List of 2
  .. ..$ x: language ~234",
      fixed = TRUE
    )
    expect_output(
      str(
        modify(
          sl1,
          X1 = s2
        )
      ),
      "List of 2
 $ X1:List of 2
  ..$ .dots          :List of 2
  .. ..$ x: language ~987",
      fixed = TRUE
    )
    
    expect_output(
      str(
        modify(
          sl1,
          X3 = s2
        )
      ),
      "List of 3
 $ X1:List of 2
  ..$ .dots          :List of 2
  .. ..$ x: language ~234",
      fixed = TRUE
    )
    expect_error(
      heemod:::define_state_list(
        X1 = s1,
        X1 = s2
      )
    )
    expect_error(
      heemod:::define_state_list(
        X1 = define_state(x = 1, y = 2),
        X2 = define_state(x = 1)
      )
    )
    expect_error(
      heemod:::define_state_list(
        X1 = define_state(x = 1, y = 2),
        X2 = define_state(x = 1, z = 2)
      )
    )
    expect_warning(
      heemod:::define_state_list(
        X1 = s1,
        s2
      )
    )
    expect_error(
      heemod:::define_state_list(
        1:2,
        s1
      )
    )
    expect_equal(
      get_state_number(sl1), 2
    )
    expect_equal(
      get_state_names(sl1),
      c("X1", "X2")
    )
    expect_equal(
      get_state_value_names(sl1),
      c("x", "y")
    )
  }
)

test_that(
  "State evaluation", {
    par1 <- define_parameters(
      a = 12,
      b = a + 4 * model_time
    )
    sl <- define_state_list(
      X1 = define_state(A = a),
      X2 = define_state(A = b)
    )
    e_par1 <- eval_parameters(
      par1, 10
    )
    e_st <- eval_state_list(sl, e_par1)
    expect_output(
      print(e_st),
      "2 evaluated states, 10 Markov cycles.",
      fixed = TRUE
    )
    expect_equal(names(e_st), c(".dots", "starting_values"))
    expect_output(
      str(e_st$.dots),
      "..$ model_time: int [1:10] 1 2 3 4 5 6 7 8 9 10
  ..$ A         :",
      fixed = TRUE
    )
    expect_equal(
      get_state_number(e_st), 2
    )
    expect_equal(
      get_state_names(e_st),
      c("X1", "X2")
    )
    expect_equal(
      get_state_value_names(e_st),
      c("A")
    )
  }
)

test_that(
  "Discounting works", {
    
    par1 <- define_parameters(
      a = .1,
      b = 1 / (model_time + 1),
      cte1 = 234,
      cte2 = 123,
      cte3 = 987,
      cte4 = 1726
    )
    mat1 <- define_transition(
      state_names = c("X1", "X2"),
      1-a, a,
      1-b, b
    )
    s1 <- define_state(
      x = discount(234, .05),
      y = discount(123, .1) * 1
    )
    f <- function(x) x
    s2 <- define_state(
      x = f(discount(987, .2)),
      y = discount(r = .1, 1726)
    )
    mod1 <- define_strategy(
      transition = mat1,
      X1 = s1,
      X2 = s2
    )
    s3 <- define_state(
      x = discount(cte1, .05),
      y = discount(cte2, .1)
    )
    s4 <- define_state(
      x = discount(cte3, .2),
      y = discount(cte4, .1)
    )
    mod2 <- define_strategy(
      transition = mat1,
      X1 = s3,
      X2 = s4
    )
    res <- run_model(
      mod1, mod2,
      parameters = par1, cost = x, effect = y
    )
    expect_equal(
      summary(res)$res_comp$.icer, c(NA, NaN)
    )
    
    discount <- function(x, r) x
    expect_warning(
      run_model(
        mod1, mod2,
        parameters = par1, cost = x, effect = y
      ),
      "version"
    ) |> 
      suppressWarnings()
  }
)

context("Environments")

f1 <- function(x) x + 1
test_that(
  "Current environment > caller environment for define_parameters", {
    age_init <- 25
    res <- eval_parameters(define_parameters(age_init = 12, age = age_init + model_time))
    expect_equal(res$age, 13)
    expect_error(eval_parameters(define_parameters(age = age_init + model_time)))
    res <- eval_parameters(define_parameters(age_init = find(age_init), age = age_init + model_time))
    expect_equal(res$age, 26)
})

test_that(
  "functions are found", {
    f2 <- function(x) x + 2
    # eval_parameters(define_parameters(age = 25,
    #                                   age_f = f(age),
    #                                   age_f2 = f2(age))) %>% 
    #   expect_error()
    params <- define_parameters(
      age = 25, 
      f1 = find(f1),
      age_f = f1(age), 
      f2 = find(f2),
      age_f2 = f2(age_f))

    res <- eval_parameters(params)
    expect_equal(res$age_f, 26)
    expect_equal(res$age_f2, 28)
  })

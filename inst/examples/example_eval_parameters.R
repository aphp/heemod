param <- define_parameters(
  age_start = 60,
  age = age_start + model_time
)

heemod:::eval_parameters(param, cycles = 15)

# heemod

# Changelog

## heemod 1.1.0
### Breaking changes
 * Remove dependency to rgho. `get_who_mr` is still locally working

## heemod 1.0.0
### Breaking changes
 * dplyr's pipe now causes an error when chaining survival operation. The new pipe |>
is required in this situation.
 * `define_survival()` has been renamed `define_surv_dist()` for naming consistency
 * It is now necessary to use `define_surv_fit()` to specify a survival fit created
 with `survfit()`, `flexsurvreg()` or `flexsurvspline()`

### New Features
 * It is now possible to run probabilistic sensitivity analysis (PSA) with survival
 objects - either created from a parametric distribution (with `define_surv_dist`)
 or with real data fitted with `survfit()` or `flexsurvreg()`.
 See `vignette("j-survival_2_psa")`
 * `plot.surv_object()` gains `psa` argument allowing to plot the distribution to
 check for its variability.
 * Created method `c.uneval_parameters()` to combine outputs of `define_parameters()`.
 
### Other changes
 * All operation functions on survival objects are now of class `surv_object`

## heemod 0.16.0
 * Remove lazyeval, plyr and pryr dependency
 * For PSM, when OS is lower to PFS, force the transition probability to be 0 (#7 @william-rawlinson)
 * `get_who_mr()` is working again
 * Mortality rates gathered with `get_who_mr(..., local = TRUE)` are only available as     
 pooled for high income countries.
 * `prop`, `multinom` and `make_gamma`, deprecated since a long time are now unavailable,
 to be replaced with `binomial`, `multinomial` and `gamma`, respectively.
 * `markov_cycle` is now deprecated in favor of `model_time` and `state_cycle` in favor of `state_time`

## heemod 0.15.1
Compliance with r-devel

## heemod 0.15.0
  * Performance improvements: using state_time is much faster now, as well as  probabilistic sensitivity analyses

## heemod 0.14.5
  * Rebuild vignettes
  * fix guides(<scale> = FALSE)

## heemod 0.14.4
  * updates for compatibility with rlang 1.0.0

## heemod 0.14.2
  * checks cairo_pdf availability

## heemod 0.14.0

### Breaking changes
  * argument "method" for run_model was inverted. Now "beginning" overestimates costs and "end" underestimates costs (#345 @zaddyzad)
  
## heemod 0.13.0

Prepare for dplyr 1.0.0

### Bugfixes
  * Fix problem with state_time (fix #343 @Reale3)

## heemod 0.12.0

### Bugfixes
  * Ready for tibble 3.0
  * `look_up()` displays a more informative message when searching several values at once

## heemod 0.11.0

### Backend changes

  * States are now length 2 lists: 
    * .dots: the state values
    * starting_values: the state values when entering the state
    
  * Compatibility with dplyr 0.8: replaced all underscored (SE) versions of all verbs
  
### Bugfixes
  
  * `get_who_mr()` is again able to fetch data from GHO server

## heemod 0.10.0

### New features

  * `define_state()` gains the `starting_values` argument (same as `define_strategy()`).

## heemod 0.9.4

  * Compatibility with dplyr v0.8

## heemod 0.9.3

### Backend changes

  * Remove `Hmisc` dependency for weighted summaries.

## heemod 0.9.2

### Breaking changes

  * Survival function renamed: `pool()`=>`mix()`.

### New features

  * Starting values with `define_starting_values()`.
  * Check for comma as a decimal separator in tabular data.
  * Auto-re-indent transition matrices with `reindent_transition()`.
  * Discounting per year with the `period` argument.
  * [BCEA](https://sites.Google.com/a/statistica.it/gianluca/bcea) interface with `run_bcea()`.
  * `use_distribution()` to re-sample from empirical data (e.g. MCMC fit posterior distribution).

### Bugfixes

  * Better error-checking in multinomial parameters definition.
  * Cleaner `define_psa()` back-end code, variables can now be used to define distribution parameters.
  * Handle cases when parameters depend on `state_time` but those parameters are not used by a strategy.
  * Empty calls to `modify()` are supported.
  * Fix error in survival vignette.

### Backend changes

  * The development branch is now called `master`, the release branch is called `cran`.
  * Survival operation are wrapped in a common class.
  * Expansion code moved to own function.
  * Faster computation of expanded states.
  * Stricter dependency version management (to avoid install failures on an outdated system).
  * Removed dependency to `tidyr`.

## heemod 0.9.1

### Breaking changes

  * Survival function renamed: `project()`=>`join()`.
  * Removed `get_code()`.

### New features

  * NMB is now reported as relative values.

### Bugfixes

  * Fixed error with `run_dsa()` in state expansion situations.

### Backend changes

  * Multiple changes to prepare for the `dplyr` update 0.6.0.

## heemod 0.9.0

`heemod` can now be cited with:

Filipovic-Pierucci A, Zarca K and Durand-Zaleski I (2017).
“Markov Models for Health Economic Evaluations: The R
Package heemod.” _ArXiv e-prints_. R package version
0.8.0, 1702.03252

  * Model calibration with `calibrate_model()`.
  * `init` and `inflow` can refer to values from `define_parameters()`, and are thus accessible to PSA / DSA.
  * `inflow` can be time-varying.

### Breaking changes

  * `get_probs_from_surv()`=>`compute_surv()`.

### New features

  * Operations on survival functions, see `vignette("j-survival")`. Thanks to [Jordan Amdahl](https://github.com/jrdnmdhl).
  * Model calibration with `calibrate_model()`. Thanks to [Matthew Wiener](https://github.com/MattWiener).
  * State membership correction works for partitioned survival model.
  * Custom state membership correction functions can be passed to the `method` argument of `run_model()`.
  * `get_who__mr()` can take a `region` argument, and automatically pool results when `sex` is not passed.
  * Additional DSA plotting options.
  
### Bugfixes

  * `heemod::discount()` now produces results similar to `discount()`.

### Backend changes

  * Renamed memoised functions.
  * NSE aliases (`*_`) to survival operation functions.

## heemod 0.8.0

### Breaking changes

  * `strategy` becomes a reserved parameter name.
  * Counting method `"half-cycle"` is deprecated.
  * `resample` argument in `run_psa()` renamed to `psa`.
  * Probability distributions for PSA were renamed, see `?distributions`.
  * `state_cycle` renamed to `state_time`, added a `model_time` alias to `model_time`.
  * The `inflow` argument for budget impact analysis is specified with `define_inflow()`.
  * `prob_to_prob()`=>`rescale_prob()`.

### New features

  * Transition probabilities from survival models with `get_probs_from_surv()`:
    * Parametrically defined by `define_survival()`.
    * Or fitted on data by the `flexsurv` package.
  * Support for partitioned survival models with `define_part_surv()`.
  * Plot EVPI.
  * Export PSA files for Sheffield Accelerated Value of Information software.
  * Individuals can enter the model after the beginning with the `inflow` argument in `run_model()` (mainly for budget impact analysis).
  * Strategy name can be used to define values with `dispatch_strategy()` or using the `strategy` name (vignettes *homogeneous* and *probabilistic* have been updated to use this feature).
  * Beta and triangle distributions for PSA.
  * Custom distributions can be defined.
  * Co-variance analysis on strategy differences, more options for `gam()` fitting.

### Other features

  * Convenience function `rescale_discount_rate()`.
  * Better error messages at parameter evaluation.
  * `combine_probs()`: given several independent probabilities of an event, return the total probability of the event.
  * More informative error messages for incorrect matrices.
  * Infinite parameter values generate an error.
  * New parameters can be added with `modify()`.

### Backend changes

  * Cleaner handling of clusters.
  * Package `diagram`, `Hmisc` and `logitnorm` moved from `Imports` to `Suggests`.
  
### Bugfixes

  * Fixed several bugs that would return incorrect efficiency frontiers in some situations, or would return duplicated strategy names in some edge cases (thanks to [Vince Daniels](https://github.com/daniels4321)).
  * Fixed a failure of tabular input when a column could be read as all numeric.
  * Character variables from `newdata` were mistakenly parsed as lazy expressions.
  * `newdata` now handles factor variables.
  * `state_cycle_limit` was not passed to PSA, DSA, or updating.
  * PSA and DSA tabular files were not saved.
  
### Acknowledgments

  * Thanks to [Matthew Wiener](https://github.com/MattWiener), especially for the survival analysis code.
  
## heemod 0.7.1

### Bugfixes

  * Fixed an error resulting in incorrect co-variance analysis results when relations between values and parameters were negative.

## heemod 0.7.0

### Breaking changes

  * Some plotting arguments changed (e.g. `type`=>`results`).

### New features

  * Added co-variance analysis for PSA with `type = "cov"` plot.
  * All plotting functions can now represent multiple strategies with facets.
  * CECA plotted on a log scale by default.
  * Black & white plots for publications with the `bw` plot option.
  * Remove variables that are not sensitive in DSA with `remove_ns`.

### Backend changes

  * Important object structure change (from attributes to lists).
  * More systematic use of getter functions.
  * More standardized processing of model objects.
  * Unit tests rely less on printed results.
  * Using new nomenclature in back-end functions (`model`=>`strategy`).
  * `base_strategy` divided in 3 concepts: `central_strategy`, `root_strategy`, `noncomparable_strategy`.
  * Non-heemod versions of `discount()` throw warnings.

## heemod 0.6.0

### Breaking changes

  * Base model cannot be specified anymore: it is always the least costly model.
  * Renamed argument `transition_matrix` => `transition` in `define_strategy()`.

### New features

  * Parallel computing with `use_cluster()`.
  * Compute average values from PSA.
  
### Bugfixes
  
  * Acceptability curve returns probabilities at 0.
  * Correctly identify efficiency frontier.
  * ICER are computed on the efficiency frontier.
  
### Acknowledgments

  * Thanks to [Matthew Wiener](https://github.com/MattWiener), [Zdenek Kabat](https://github.com/zkabat) and [Vojtech Filipec](https://github.com/vojtech-filipec) for their great contributions to this update.

## heemod 0.5.1

### Bugfixes

  * Added `formatR` to suggest list.

## heemod 0.5.0

### Breaking changes

  * Some functions were renamed for clarification:
    *  `define_model()` => `define_strategy()`
    *  `run_models()` => `run_model()`
    *  `define_matrix()` => `define_transition()`
    *  `define_sensitivity()` => `define_dsa()`
    *  `define_distrib()` => `define_psa()`
    *  `run_sensitivity()` => `run_dsa()`
    *  `run_probabilistic()` => `run_psa()`

### New features

  * Values and probabilities can depend on state time with `state_cycle`, allowing to reproduce the results of micro-simulations.
  * `define_sensitivity()` now accepts any expression as input, and can call references to model parameters.
  * Discount rates can now be specified as parameters (allows for rates to be modified in DSA & PSA).
  * Any state value can be plotted.
  * Additional output: csv files for many of the tabular results.
  
### Bugfixes

  * Edges where P=0 are not plotted anymore for transition matrices.
  * There used to be a potential error when you wanted to save output but there had been no psa.

## heemod 0.4.0

### New features

  * `update()` for heterogeneity analysis and to compute population-level values, with vignette.
  * `run_models_tabular()` to import models from tabular input, with vignette.
  * `look_up()` to look up values from external data.
  * Added option to pool female and male mortality rates in WHO data.
  * Counting method now defaults to life-table.
  * `plot_sensitivity()` now plots by default the widest bar on top.
  * Convenience functions for converting rates to probabilities.
  * Models can be run without state values, to compute counts only.
  * Much more informative error messages.
  * Objects can be converted to the `R` code to generate them (same idea as `dput()`, but easier to read).
  * New options `heemod.verbose` and `heemod.memotime`.
  * More informative messages, especially in verbose mode.
  * Use WHO data cached locally in case of connection problems.
  * New functions: `get_counts()` and `get_init()` to get state membership counts.
  * Smart sex code conversion for `get_who_mr()`.

### Bug fixes

  * _really_ fixed problem when the argument to `discount()` was not defined as a parameter.
  
### Backend changes

  * `eval_model_newdata()`, the function behind re-sampling and sensitivity analysis now returns list-variables.
  
### Acknowledgments

  * Thanks to [Matthew Wiener](https://github.com/MattWiener), [Zdenek Kabat](https://github.com/zkabat) and [Vojtech Filipec](https://github.com/vojtech-filipec) for their great contributions to this update.

## heemod 0.3.3

### New features

  * Added a vignette to exactly reproduce results from Decision Modelling for Health Economic Evaluation.

### Bug fixes

  * Fixed problem when argument to `discount()` was not defined as a parameter.
  * Corrected several errors in the vignettes (thanks to [Michael Schenkenberg](https://github.com/MichaelSvm) from SBU, Stockholm, Sweden).
  * Updated mortality rate tests to reflect GHO database update.

## heemod 0.3.2

### Bug fixes

  * Fix mishandling of matrix index with `C` in `eval_matrix()` (thanks to [Matthew Wiener](https://github.com/MattWiener)).
  * Fix problem with upcoming version of `tidyr`.

## heemod 0.3.1

### Backend changes

  * A single set of parameters is now used for a set of models.

### Bug fixes

  * Correct error in probabilistic analysis vignette.
  * Fix incompatibility with upcoming `dplyr` update.

## heemod 0.3.0

### New features

  * `shiny` interface.
  * Added support for different counting options.
  * `get_who_mr()` to extract mortality rates from WHO data.

## heemod 0.2.0

### New features

  * Added ability to plot matrix.
  * Added framework for probabilistic uncertainty analysis.
  * Written vignette for probabilistic and sensitivity analysis.
  * Added complement alias for matrix definition.
  * `run_*` functions now output a single table instead of a list of tables.
  * Variables corresponding to cost and effect must now be specified in `run_models()`.
  * Plotting for sensitivity and probabilistic analysis.
  
### Bug fixes

  * added `*_` functions.
  
### Removed

  * `run_model_newdata()` cannot run `uneval_model` anymore.
  * Removed need to use `define_state_list()`.

## heemod 0.1.0

  * Initial CRAN submission.

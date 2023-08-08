## Evaluation of the variability of the survival distribution
\donttest{
surv1 <- define_surv_dist("exp", rate = 0.1)
psa <- define_psa(surv1 ~ resample_surv(n = 100))
plot(surv1, psa=psa)

## plot surv_projection object
surv2 <- define_surv_dist("exp", rate = 0.5)
plot(join(surv1, surv2, at = 2), psa = psa, Nrep = 50)

## surv_fit object
library(survival)
km <- define_surv_fit(survfit(formula = Surv(time, status) ~ 1, data = aml))
fs <- flexsurv::flexsurvreg(formula = Surv(time, status) ~ 1, 
                        data = aml, 
                        dist = "weibull") |>
  define_surv_fit()

psa2 <- define_psa(km ~ resample_surv(),
                   fs ~ resample_surv(),
                   surv1 ~ resample_surv(100))
plot(km, psa = psa2)

plot(join(km, surv1, at = 6), psa = psa2)
plot(join(fs, surv1, at = 6), psa = psa2)
}

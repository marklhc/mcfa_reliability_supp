# Data Import and Setup ---------------------------------------------------

library(tidyverse)

timss15_eng <- here::here("timss/data", "MSGUSAM3.sav") %>% 
  haven::read_sav() %>% 
  select(IDCNTRY, IDSCHOOL, IDCLASS, 
         MSBM18A, MSBM18B, MSBM18C, MSBM18H, MSBM18I, MSBM18M)

# Export to Mplus
timss15_eng %>% 
  select(IDCLASS, 
         MSBM18A, MSBM18B, MSBM18C, MSBM18H, MSBM18I, MSBM18M) %>% 
  write_delim(here::here("timss", "timss15_eng.dat"), 
              na = "-99", col_names = FALSE)

# Harmonic mean of school sizes (25.1):
count(timss15_eng %>% drop_na(), IDCLASS) %>% 
  transmute(n_inv = 1 / n) %>% 
  summarize(1 / mean(n_inv))

# Obtain omega^w and 95% CI using a saturated between-cluster model --------
# For a within-cluster construct

library(lavaan)
mcfas1 <- 
  'level: 1
     MSBM18A ~~ c11w * MSBM18A + c12w * MSBM18B + c13w * MSBM18C + 
                c14w * MSBM18H + c15w * MSBM18I + c16w * MSBM18M
     MSBM18B ~~ c22w * MSBM18B + c23w * MSBM18C + c24w * MSBM18H + 
                c25w * MSBM18I + c26w * MSBM18M
     MSBM18C ~~ c33w * MSBM18C + c34w * MSBM18H + c35w * MSBM18I + 
                c36w * MSBM18M
     MSBM18H ~~ c44w * MSBM18H + c45w * MSBM18I + c46w * MSBM18M
     MSBM18I ~~ c55w * MSBM18I + c56w * MSBM18M
     MSBM18M ~~ c66w * MSBM18M
   level: 2
     f1b =~ NA * MSBM18A + l1 * MSBM18A + l2 * MSBM18B + l3 * MSBM18C + 
            l4 * MSBM18H + l5 * MSBM18I + l6 * MSBM18M
     # fixed residual variances
     MSBM18A ~~ ev1b * MSBM18A
     MSBM18B ~~ ev2b * MSBM18B
     MSBM18C ~~ ev3b * MSBM18C
     MSBM18H ~~ ev4b * MSBM18H
     MSBM18I ~~ ev5b * MSBM18I
     MSBM18M ~~ ev6b * MSBM18M
     f1b ~~ 1 * f1b
   # score reliability:
   omgb := (l1 + l2 + l3 + l4 + l5 + l6)^2 /
              ((l1 + l2 + l3 + l4 + l5 + l6)^2 + 
               ev1b + ev2b + ev3b + ev4b + ev5b + ev6b + 
           (c11w + c22w + c33w + c44w + c55w + c66w + 
            2 * (c12w + c13w + c14w + c15w + c16w + 
                 c23w + c24w + c25w + c26w + 
                 c34w + c35w + c36w + c45w + c46w + c56w)) / 9.15)'

mcfas1_fit <- cfa(mcfas1, data = timss15_eng, cluster = "IDCLASS")
summary(mcfas1_fit, fit.measures = TRUE, ci = TRUE)

# Obtain Monte Carlo CI:
coef_to_extract <- c("l1", "l2", "l3", "l4", "l5", "l6", 
                     "ev1b", "ev2b", "ev3b", "ev4b", "ev5b", "ev6b", 
                     "c11w", "c22w", "c33w", "c44w", "c55w", "c66w", 
                     "c12w", "c13w", "c14w", "c15w", "c16w", 
                     "c23w", "c24w", "c25w", "c26w", 
                     "c34w", "c35w", "c36w", "c45w", "c46w", "c56w")
est_mcfas1 <- coef(mcfas1_fit)[coef_to_extract]
acov_mcfas1 <- vcov(mcfas1_fit)[coef_to_extract, coef_to_extract]
# Simulate new coefficients based on a multivariate normal distribution
ndraws <- 1e4
set.seed(123)
simcoef_mcfas1 <- MASS::mvrnorm(ndraws, mu = est_mcfas1, 
                                Sigma = acov_mcfas1)
# Compute omegab
sim_omgb <- with(data.frame(simcoef_mcfas1), 
                 (l1 + l2 + l3 + l4 + l5 + l6)^2 /
                   ((l1 + l2 + l3 + l4 + l5 + l6)^2 + 
                      ev1b + ev2b + ev3b + ev4b + ev5b + ev6b + 
                      (c11w + c22w + c33w + c44w + c55w + c66w + 
                         2 * (c12w + c13w + c14w + c15w + c16w + 
                                c23w + c24w + c25w + c26w + 
                                c34w + c35w + c36w + c45w + c46w + c56w)) / 9.15))
# Histogram
hist(sim_omgb)
# 95% Monte Carlo CI
quantile(sim_omgb, probs = c(.025, .975))

# Simultaneous shared-and-configural model --------------------------------

# Run MCFA (configural + shared):
mcfa11shared <- 
  'level: 1
     f1w =~ NA * MSBM18A + l1 * MSBM18A + l2 * MSBM18B + l3 * MSBM18C + 
            l4 * MSBM18H + l5 * MSBM18I + l6 * MSBM18M
     MSBM18A ~~ ev1w * MSBM18A
     MSBM18B ~~ ev2w * MSBM18B
     MSBM18C ~~ ev3w * MSBM18C
     MSBM18H ~~ ev4w * MSBM18H
     MSBM18I ~~ ev5w * MSBM18I
     MSBM18M ~~ ev6w * MSBM18M
     f1w ~~ 1 * f1w
   level: 2
     f1b =~ NA * MSBM18A + l1 * MSBM18A + l2 * MSBM18B + l3 * MSBM18C + 
            l4 * MSBM18H + l5 * MSBM18I + l6 * MSBM18M
     f1bs =~ MSBM18A + l2s * MSBM18B + l3s * MSBM18C + 
             l4s * MSBM18H + l5s * MSBM18I + l6s * MSBM18M
     # fixed residual variances
     MSBM18A ~~ ev1b * MSBM18A
     MSBM18B ~~ ev2b * MSBM18B
     MSBM18C ~~ ev3b * MSBM18C
     MSBM18H ~~ ev4b * MSBM18H
     MSBM18I ~~ ev5b * MSBM18I
     MSBM18M ~~ ev6b * MSBM18M
     f1b ~~ (.05 / (1 - .05)) * f1b
     f1bs ~~ v1bs * f1bs
     f1b ~~ 0 * f1bs
   # score reliability:
   omgbs := (1 + l2s + l3s + l4s + l5s + l6s)^2 * v1bs /
            ((l1 + l2 + l3 + l4 + l5 + l6)^2 * .05 / (1 - .05) + 
             (1 + l2s + l3s + l4s + l5s + l6s)^2 * v1bs + 
             ev1b + ev2b + ev3b + ev4b + ev5b + ev6b + 
             (ev1w + ev2w + ev3w + ev4w + ev5w + ev6w + 
              (l1 + l2 + l3 + l4 + l5 + l6)^2) / 9.15)'
mcfa11sh_fit <- cfa(mcfa11shared, data = timss15_eng, cluster = "IDCLASS")
summary(mcfa11sh_fit, fit.measures = TRUE, ci = TRUE)

# Obtain Monte Carlo CI:
coef_to_extract <- c("l1", "l2", "l3", "l4", "l5", "l6", 
                     "l2s", "l3s", "l4s", "l5s", "l6s", 
                     "ev1w", "ev2w", "ev3w", "ev4w", "ev5w", "ev6w",
                     "ev1b", "ev2b", "ev3b", "ev4b", "ev5b", "ev6b",
                     "v1bs")
est_mcfa11sh <- coef(mcfa11sh_fit)[coef_to_extract]
acov_mcfa11sh <- vcov(mcfa11sh_fit)[coef_to_extract, coef_to_extract]
# Simulate new coefficients based on a multivariate normal distribution
ndraws <- 1e4
set.seed(123)
simcoef_mcfa11sh <- MASS::mvrnorm(ndraws, mu = est_mcfa11sh, 
                                  Sigma = acov_mcfa11sh)
# Compute omega_bs (shared component at between level)
sim_omgbs <- with(data.frame(simcoef_mcfa11sh), 
                  (1 + l2s + l3s + l4s + l5s + l6s)^2 * v1bs /
                    ((l1 + l2 + l3 + l4 + l5 + l6)^2 * .05 / (1 - .05) + 
                       (1 + l2s + l3s + l4s + l5s + l6s)^2 * v1bs + 
                       ev1b + ev2b + ev3b + ev4b + ev5b + ev6b + 
                       (ev1w + ev2w + ev3w + ev4w + ev5w + ev6w + 
                          (l1 + l2 + l3 + l4 + l5 + l6)^2) / 9.15))
# Histogram
hist(sim_omgbs)
# 95% Monte Carlo CI
quantile(sim_omgbs, probs = c(.025, .975))

# Data Import and Setup ---------------------------------------------------

library(tidyverse)

timss_usa <- here::here("timss/data", "asgUSAm4.sav") %>% 
   haven::read_sav() %>% 
   select(IDCNTRY, IDSCHOOL, 
          AS4MAMOR, AS4MAENJ, AS4MALIK, AS4MABOR) %>% 
   drop_na() %>% 
   mutate(AS4MABORr = 5 - AS4MABOR)

# Harmonic mean of school sizes (25.1):
count(timss_usa, IDSCHOOL) %>% 
   transmute(n_inv = 1 / n) %>% 
   summarize(1 / mean(n_inv))

# Run MCFA (assume continuous items; for an individual construct):
library(lavaan)
mcfa11 <- 
   'level: 1
     f1w =~ NA * AS4MAMOR + l1 * AS4MAMOR + l2 * AS4MAENJ + l3 * AS4MALIK + l4 * AS4MABORr
     AS4MAMOR ~~ ev1w * AS4MAMOR
     AS4MAENJ ~~ ev2w * AS4MAENJ
     AS4MALIK ~~ ev3w * AS4MALIK
     AS4MABORr ~~ ev4w * AS4MABORr
     f1w ~~ 1 * f1w
   level: 2
     f1b =~ NA * AS4MAMOR + l1 * AS4MAMOR + l2 * AS4MAENJ + l3 * AS4MALIK + l4 * AS4MABORr
     # fixed residual variances
     AS4MAMOR ~~ ev1b * AS4MAMOR
     AS4MAENJ ~~ 0 * AS4MAENJ
     AS4MALIK ~~ ev3b * AS4MALIK
     AS4MABORr ~~ ev4b * AS4MABORr
     f1b ~~ vf1b * f1b
   # tilde omega values:
   tilomgb := (l1 + l2 + l3 + l4)^2 * vf1b /
              ((l1 + l2 + l3 + l4)^2 * vf1b + ev1b + 0 + ev3b + ev4b)
   tilomgw := (l1 + l2 + l3 + l4)^2 * 1 /
              ((l1 + l2 + l3 + l4)^2 * 1 + ev1w + ev2w + ev3w + ev4w)
   # score reliability:
   omg2l := (l1 + l2 + l3 + l4)^2 * (1 + vf1b) /
            ((l1 + l2 + l3 + l4)^2 * (1 + vf1b) + 
             ev1b + 0 + ev3b + ev4b + ev1w + ev2w + ev3w + ev4w)
   omgb := (l1 + l2 + l3 + l4)^2 * vf1b /
           ((l1 + l2 + l3 + l4)^2 * vf1b + ev1b + 0 + ev3b + ev4b + 
            (ev1w + ev2w + ev3w + ev4w + (l1 + l2 + l3 + l4)^2) / 25.1)'

mcfa11_fit <- cfa(mcfa11, data = timss_usa, cluster = "IDSCHOOL")
summary(mcfa11_fit, fit.measures = TRUE)

# Obtain Monte Carlo CI:
coef_to_extract <- c("l1", "l2", "l3", "l4", 
                     "ev1w", "ev2w", "ev3w", "ev4w", 
                     "ev1b", "ev3b", "ev4b", 
                     "vf1b")
est_mcfa11 <- coef(mcfa11_fit)[coef_to_extract]
acov_mcfa11 <- vcov(mcfa11_fit)[coef_to_extract, coef_to_extract]
# Simulate new coefficients based on a multivariate normal distribution
ndraws <- 1e4
set.seed(123)
simcoef_mcfa11 <- MASS::mvrnorm(ndraws, mu = est_mcfa11, 
                                Sigma = acov_mcfa11)
# Compute omega for each simulation
sim_omg2l <- with(data.frame(simcoef_mcfa11), 
                  (l1 + l2 + l3 + l4)^2 * (1 + vf1b) /
                     ((l1 + l2 + l3 + l4)^2 * (1 + vf1b) + 
                         ev1b + 0 + ev3b + ev4b + ev1w + ev2w + ev3w + ev4w))
# Compute omega_w for each simulation
sim_omgw <- with(data.frame(simcoef_mcfa11), 
                 (l1 + l2 + l3 + l4)^2 * 1 /
                    ((l1 + l2 + l3 + l4)^2 * 1 + ev1w + ev2w + ev3w + ev4w))
# Compute omega_b for each simulation
sim_omgb <- with(data.frame(simcoef_mcfa11), 
                 (l1 + l2 + l3 + l4)^2 * vf1b /
                    ((l1 + l2 + l3 + l4)^2 * vf1b + ev1b + 0 + ev3b + ev4b + 
                        (ev1w + ev2w + ev3w + ev4w + 
                            (l1 + l2 + l3 + l4)^2) / 25.1))
# Histogram
hist(sim_omg2l)
hist(sim_omgw)
hist(sim_omgb)
# 95% Monte Carlo CI for omega^2l and omega^w
quantile(sim_omg2l, probs = c(.025, .975))
quantile(sim_omgw, probs = c(.025, .975))
quantile(sim_omgb, probs = c(.025, .975))


# Obtain omega^w and 95% CI using a saturated between-cluster model --------
# For a within-cluster construct

mcfa1s <- 
   'level: 1
     f1w =~ NA * AS4MAMOR + l1 * AS4MAMOR + l2 * AS4MAENJ + l3 * AS4MALIK + l4 * AS4MABORr
     # fixed residual variances
     AS4MAMOR ~~ ev1w * AS4MAMOR
     AS4MAENJ ~~ ev2w * AS4MAENJ
     AS4MALIK ~~ ev3w * AS4MALIK
     AS4MABORr ~~ ev4w * AS4MABORr
     f1w ~~ 1 * f1w
    level: 2
     AS4MAMOR ~~ c11b * AS4MAMOR + c12b * AS4MAENJ + 
                 c13b * AS4MALIK + c14b * AS4MABORr
     AS4MAENJ ~~ c22b * AS4MAENJ + c23b * AS4MALIK + c24b * AS4MABORr
     AS4MALIK ~~ c33b * AS4MALIK + c34b * AS4MABORr
     AS4MABORr ~~ c44b * AS4MABORr
   # tilde omega values:
   omgw := (l1 + l2 + l3 + l4)^2 /
              ((l1 + l2 + l3 + l4)^2 + ev1w + ev2w + ev3w + ev4w)'

mcfa1s_fit <- cfa(mcfa1s, data = timss_usa, cluster = "IDSCHOOL")
summary(mcfa1s_fit, fit.measures = TRUE)

# Obtain Monte Carlo CI:
coef_to_extract <- c("l1", "l2", "l3", "l4",  
                     "ev1w", "ev2w", "ev3w", "ev4w")
est_mcfa1s <- coef(mcfa1s_fit)[coef_to_extract]
acov_mcfa1s <- vcov(mcfa1s_fit)[coef_to_extract, coef_to_extract]
# Simulate new coefficients based on a multivariate normal distribution
ndraws <- 1e4
set.seed(123)
simcoef_mcfa1s <- MASS::mvrnorm(ndraws, mu = est_mcfa1s, 
                                Sigma = acov_mcfa1s)
# Compute omegab
sim_omgw <- with(data.frame(simcoef_mcfa1s), 
                 (l1 + l2 + l3 + l4)^2 /
                    ((l1 + l2 + l3 + l4)^2 + ev1w + ev2w + ev3w + ev4w))
# Histogram
hist(sim_omgw)
# 95% Monte Carlo CI
quantile(sim_omgw, probs = c(.025, .975))

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

# Run MCFA (assume continuous items):
library(lavaan)
msat <- 
  'level: 1
     AS4MAMOR ~~ t11w * AS4MAMOR + t12w * AS4MAENJ + t13w * AS4MALIK + t14w * AS4MABORr
     AS4MAENJ ~~ t22w * AS4MAENJ + t23w * AS4MALIK + t24w * AS4MABORr
     AS4MALIK ~~ t33w * AS4MALIK + t34w * AS4MABORr
     AS4MABORr ~~ t44w * AS4MABORr
   level: 2
     AS4MAMOR ~~ t11b * AS4MAMOR + t12b * AS4MAENJ + t13b * AS4MALIK + t14b * AS4MABORr
     AS4MAENJ ~~ t22b * AS4MAENJ + t23b * AS4MALIK + t24b * AS4MABORr
     AS4MALIK ~~ t33b * AS4MALIK + t34b * AS4MABORr
     AS4MABORr ~~ t44b * AS4MABORr
   # tilde alpha values:
   tilalpb := 4 * 2 * (t12b + t13b + t14b + t23b + t24b + t34b) / 
              (4 - 1) / (t11b + t22b + t33b + t44b + 
                         2 * (t12b + t13b + t14b + t23b + t24b + t34b))
   tilalpw := 4 * 2 * (t12w + t13w + t14w + t23w + t24w + t34w) / 
              (4 - 1) / (t11w + t22w + t33w + t44w + 
                         2 * (t12w + t13w + t14w + t23w + t24w + t34w))
   # score reliability:
   alpha2l := 4 * 2 * (t12b + t13b + t14b + t23b + t24b + t34b + 
                       t12w + t13w + t14w + t23w + t24w + t34w) / 
              (4 - 1) / (t11b + t22b + t33b + t44b + 
                         2 * (t12b + t13b + t14b + t23b + t24b + t34b) + 
                         t11w + t22w + t33w + t44w + 
                         2 * (t12w + t13w + t14w + t23w + t24w + t34w))
   alphab := 4 * 2 * (t12b + t13b + t14b + t23b + t24b + t34b) / 
             (4 - 1) / (t11b + t22b + t33b + t44b + 
                        2 * (t12b + t13b + t14b + t23b + t24b + t34b) + 
                        (t11w + t22w + t33w + t44w + 
                         2 * (t12w + t13w + t14w + t23w + t24w + t34w)) / 25.1)
  '

msat_fit <- cfa(msat, data = timss_usa, cluster = "IDSCHOOL")
summary(msat_fit, fit.measures = TRUE)
# Obtain Monte Carlo CI:
coef_to_extract <- c("t11w", "t22w", "t33w", "t44w", 
                     "t12w", "t13w", "t14w",
                     "t23w", "t24w", 
                     "t34w", 
                     "t11b", "t22b", "t33b", "t44b", 
                     "t12b", "t13b", "t14b",
                     "t23b", "t24b", 
                     "t34b")
est_msat <- coef(msat_fit)[coef_to_extract]
acov_msat <- vcov(msat_fit)[coef_to_extract, coef_to_extract]
# Simulate new coefficients based on a multivariate normal distribution
ndraws <- 1e4
simcoef_msat <- MASS::mvrnorm(ndraws, mu = est_msat, 
                              Sigma = acov_msat)
# Compute alpha from simulated data
sim_alp2l <- with(data.frame(simcoef_msat), 
                  (4 * 2 * (t12b + t13b + t14b + t23b + t24b + t34b + 
                               t12w + t13w + t14w + t23w + t24w + t34w) / 
                      (4 - 1) / (t11b + t22b + t33b + t44b + 
                                    2 * (t12b + t13b + t14b + 
                                            t23b + t24b + t34b) + 
                                    t11w + t22w + t33w + t44w + 
                                    2 * (t12w + t13w + t14w + 
                                            t23w + t24w + t34w))))
# Compute alpha_b from simulated data
sim_alpb <- with(data.frame(simcoef_msat), 
                 (4 * 2 * (t12b + t13b + t14b + t23b + t24b + t34b) / 
                    (4 - 1) / (t11b + t22b + t33b + t44b + 
                                 2 * (t12b + t13b + t14b + t23b + t24b + t34b) + 
                                 (t11w + t22w + t33w + t44w + 
                                    2 * (t12w + t13w + t14w + 
                                           t23w + t24w + t34w)) / 25.1)))
# Compute alpha_w from simulated data
sim_alpw <- with(data.frame(simcoef_msat), 
                 4 * 2 * (t12w + t13w + t14w + t23w + t24w + t34w) / 
                    (4 - 1) / (t11w + t22w + t33w + t44w + 
                                  2 * (t12w + t13w + t14w + t23w + t24w + t34w)))
# Histogram
hist(sim_alp2l)
hist(sim_alpb)
hist(sim_alpw)
# 95% Monte Carlo CI for alpha^2l, alpha^b, and alpha^w
quantile(sim_alp2l, probs = c(.025, .975))
quantile(sim_alpb, probs = c(.025, .975))
quantile(sim_alpw, probs = c(.025, .975))
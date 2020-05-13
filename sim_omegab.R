# Generate data with very low reliability, but the theoretical "between" 
# reliability is high

library(lavaan)  # load the `lavaan` package


# Set up Simulation Parameters --------------------------------------------

J <- 10000  # number of clusters
CS <- 10  # cluster size
ICC <- .2  # intraclass correlation of latent variable
NP <- 5  # number of items
LAMBDA <- rep(.5, NP)  # factor loadings (equal across items)
THETA_B <- 0.1  # between-cluster uniqueness
THETA_W <- 1  # within-level uniqueness

# Syntax for MCFA in lavaan, with cross-level loading invariance
# The latent variable at the within-level is set to have variance 1.0
MCFA_MODEL <- 
  ' level: 1
      fw =~ NA * y1 + l1 * y1 + l2 * y2 + l3 * y3 + l4 * y4 + l5 * y5
      fw ~~ 1 * fw
      y1 ~~ ey1 * y1
      y2 ~~ ey2 * y2
      y3 ~~ ey3 * y3
      y4 ~~ ey4 * y4
      y5 ~~ ey5 * y5
    level: 2
      fb =~ NA * y1 + l1 * y1 + l2 * y2 + l3 * y3 + l4 * y4 + l5 * y5
      fb ~~ vb_t * fb
      y1 ~~ evb1 * y1
      y2 ~~ evb2 * y2
      y3 ~~ evb3 * y3
      y4 ~~ evb4 * y4
      y5 ~~ evb5 * y5
    # level-2 reliability
      omgb := (l1 + l2 + l3 + l4 + l5)^2 * vb_t /
              ((l1 + l2 + l3 + l4 + l5)^2 * vb_t + 
               evb1 + evb2 + evb3 + evb4 + evb5)
  '

# Generate Data -----------------------------------------------------------

# Function to standardize the generated latent variables
stdize <- function(x) {
  len_x <- length(x)
  (x - mean(x)) / (len_x - 1) * len_x
}

gid <- rep(seq_len(J), each = CS)  # cluster ID
set.seed(1733)  # set seed to make results reproducible
eta_b <- stdize(rnorm(J)) * sqrt(ICC / (1 - ICC))  # between-cluster latent factor
eta_w <- stdize(rnorm(J * CS))  # within-level latent factor
# Between-cluster unique factors
eb <- MASS::mvrnorm(J, mu = rep(0, NP), Sigma = diag(THETA_B, NP), 
                    empirical = TRUE)  # force means = 0 and variance = Sigma
# Between-cluster latent cluster means
ylb <- tcrossprod(eta_b, LAMBDA) + eb
# Within-level unique factors
ew <- MASS::mvrnorm(J * CS, mu = rep(0, NP), Sigma = diag(THETA_W, NP), 
                    empirical = TRUE)  # force means = 0 and variance = Sigma
# Within-level deviation from cluster means
ylw <- tcrossprod(eta_w, LAMBDA) + ew
# Combine ylb and ylw to form observed y
y <- ylb[gid, ] + ylw
colnames(y) <- paste0("y", 1:NP)
# Make y and gid a data frame
df <- data.frame(y, gid)


# Compute Reliabilities ---------------------------------------------------

# True cluster-level reliability of observed scores:
Zb <- tapply(rowSums(y), gid, mean)  # compute composite scores of cluster means
vb_true <- var(eta_b * sum(LAMBDA))  # between-cluster true score variance
(relb_true <- vb_true / var(Zb))
# Same as squared correlation between Zb and eta_b
(cor(Zb, eta_b)^2)

# Tilde-omega^b estimate from running lavaan
fit <- cfa(MCFA_MODEL, data = df, cluster = "gid")
(omgb_est <- coef(fit, type = "user")["omgb"])
# So it confirms that tilde-omega^b estimates the squared correlation between
# Zlb and the latent composite
Zlb <- rowSums(ylb)
(cor(Zlb, eta_b)^2)

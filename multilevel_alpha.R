#  R function multilevel_alpha(), a wrapper for computing 
#  the reliability indices discussed in 
#  Lai, M. H. C. (2021). Composite reliability of multilevel data: 
#      It’s about observed scores and construct meanings. 
#      Psychological Methods, 26(1), 90–102. 
#      https://doi.org/10.1037/met0000287
#  Copyright (C) 2021 Lai, Hok Chio (Mark)
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

multilevel_alpha <- function(data, id, nsim = 5000, conf_level = .95, 
                             se = "robust.huber.white", 
                             pa_args = list(fa = "pc"), ...) {
   if (!require(lavaan)) stop("The lavaan package needs to be installed.")
   if (!require(psych)) stop("The psych package needs to be installed.")
   nitem <- ncol(data)
   ynames <- paste0("y", seq_len(nitem))
   colnames(data) <- ynames
   data <- cbind(data, id = id)
   tab_id <- table(id)
   hmean_cs <- 1 / mean(1 / tab_id[tab_id > 0])
   # Alpha
   # Generate syntax for saturated model
   sat_syntax <- (function(y) {
      if (length(y) <= 1) {
         return(NULL)
      }
      paste(
         c(paste(y[1], "~~", paste(y[-1], collapse = " + ")), 
           Recall(y[-1])), 
         collapse = "\n  "
      )
   })(ynames)
   msat <- paste0("level: 1\n  ", sat_syntax, "\nlevel: 2\n  ", sat_syntax)
   msat_fit <- cfa(msat, data = data, cluster = "id", se = se, 
                   test = "none", h1 = FALSE, baseline = FALSE, ...)
   coef_msat <- coef(msat_fit, type = "user")
   vcov_msat <- vcov(msat_fit)
   vw <- names(coef_msat)[
      with(msat_fit@ParTable, which(op == "~~" & lhs == rhs & level == 1))]
   cvw <- names(coef_msat)[
      with(msat_fit@ParTable, which(op == "~~" & lhs != rhs & level == 1))]
   vb <- names(coef_msat)[
      with(msat_fit@ParTable, which(op == "~~" & lhs == rhs & level == 2))]
   cvb <- names(coef_msat)[
      with(msat_fit@ParTable, which(op == "~~" & lhs != rhs & level == 2))]
   Sw <- sum(coef_msat[vw], 2 * coef_msat[cvw])
   Sb <- sum(coef_msat[vb], 2 * coef_msat[cvb])
   alpha_const <- nitem / (nitem - 1)
   alphaw <- alpha_const * sum(2 * coef_msat[cvw]) / Sw
   alpha2l <- alpha_const * sum(2 * coef_msat[c(cvw, cvb)]) / (Sb + Sw)
   alphab <- alpha_const * sum(2 * coef_msat[cvb]) / (Sb + Sw / hmean_cs)
   sim_coef_msat <- MASS::mvrnorm(nsim, 
                                  mu = coef_msat[c(vw, vb, cvw, cvb)], 
                                  Sigma = vcov_msat[c(vw, vb, cvw, cvb), 
                                                    c(vw, vb, cvw, cvb)])
   sim_Sw <- rowSums(cbind(sim_coef_msat[ , vw], 2 * sim_coef_msat[ , cvw]))
   sim_Sb <- rowSums(cbind(sim_coef_msat[ , vb], 2 * sim_coef_msat[ , cvb]))
   sim_alphaw <- alpha_const * rowSums(2 * sim_coef_msat[ , cvw]) / sim_Sw
   sim_alpha2l <- alpha_const * rowSums(2 * sim_coef_msat[ , c(cvw, cvb)]) / 
      (sim_Sb + sim_Sw)
   sim_alphab <- alpha_const * rowSums(2 * sim_coef_msat[ , cvb]) / 
      (sim_Sb + sim_Sw / hmean_cs)
   sim_alpha_cis <- lapply(list(alpha2l = sim_alpha2l, 
                                alphab = sim_alphab, 
                                alphaw = sim_alphaw), 
                           quantile, 
                           probs = .5 + c(- conf_level, conf_level) / 2)
   # Omega
   loading_labels <- paste0("l", seq_len(nitem))
   g_syntax <- paste(loading_labels, "*", ynames, 
                     collapse = " + ")
   mcfa <- paste0("level: 1\n  fw =~ ", 
                  g_syntax, 
                  "\nlevel: 2\n  fb =~ ", 
                  g_syntax)
   mcfa_fit <- cfa(mcfa, data = data, cluster = "id", se = se, 
                   test = "none", h1 = TRUE, baseline = FALSE, ...)
   mcfa_pt <- partable(mcfa_fit)
   coef_mcfa <- coef(mcfa_fit)
   vcov_mcfa <- vcov(mcfa_fit)
   coef_mcfa[loading_labels]
   ld <- names(coef_mcfa)[with(mcfa_pt, free[which(op == "=~" &
                                                      level == 1)])]
   evw <- names(coef_mcfa)[with(mcfa_pt,
                                free[which(op == "~~" &
                                              lhs == rhs &
                                              lhs != "fw" &
                                              level == 1)])]
   fvw <- names(coef_mcfa)[with(mcfa_pt, free[which(op == "~~" &
                                                       lhs == "fw")])]
   evb <- names(coef_mcfa)[with(mcfa_pt,
                                free[which(op == "~~" &
                                              lhs == rhs & lhs != "fb" & 
                                              level == 2)])]
   fvb <- names(coef_mcfa)[with(mcfa_pt, free[which(op == "~~" &
                                                       lhs == "fb")])]
   sumldsq <- sum(1, coef_mcfa[ld])^2
   sumevw <- sum(coef_mcfa[evw])
   sumevb <- sum(coef_mcfa[evb])
   omegaw <- sumldsq * coef_mcfa[[fvw]] / 
      (sumldsq * coef_mcfa[[fvw]] + sumevw)
   omega2l <- sum(sumldsq * coef_mcfa[c(fvw, fvb)]) / 
      sum(sumldsq * coef_mcfa[c(fvw, fvb)], sumevw, sumevb)
   omegab <- sumldsq * coef_mcfa[[fvb]] / 
      (sumldsq * (coef_mcfa[[fvb]] + coef_mcfa[[fvw]] / hmean_cs) + 
          sumevb + sumevw / hmean_cs)
   sim_coef_mcfa <- MASS::mvrnorm(nsim, 
                                  mu = coef_mcfa[c(ld, fvw, fvb, evw, evb)], 
                                  Sigma = vcov_mcfa[c(ld, fvw, fvb, evw, evb), 
                                                    c(ld, fvw, fvb, evw, evb)])
   sim_sumldsq <- (1 + rowSums(sim_coef_mcfa[ , ld]))^2
   sim_sumevw <- rowSums(sim_coef_mcfa[ , evw])
   sim_sumevb <- rowSums(sim_coef_mcfa[ , evb])
   sim_omegaw <- sim_sumldsq * sim_coef_mcfa[ , fvw] / 
      (sim_sumldsq * sim_coef_mcfa[ , fvw] + sim_sumevw)
   sim_omega2l <- rowSums(sim_sumldsq * sim_coef_mcfa[ , c(fvw, fvb)]) / 
      (rowSums(sim_sumldsq * sim_coef_mcfa[ , c(fvw, fvb)]) + 
          sim_sumevw + sim_sumevb)
   sim_omegab <- sim_sumldsq * sim_coef_mcfa[ , fvb] / 
      (sim_sumldsq * (sim_coef_mcfa[ , fvb] + 
                         sim_coef_mcfa[ , fvw] / hmean_cs) + 
          sim_sumevb + sim_sumevw / hmean_cs)
   sim_omega_cis <- lapply(list(omega2l = sim_omega2l, 
                                omegab = sim_omegab, 
                                omegaw = sim_omegaw), 
                           quantile, 
                           probs = .5 + c(- conf_level, conf_level) / 2)
   # resid_corb <- resid(mcfa_fit, type = "cor")$id$cov
   # diag(resid_corb) <- 1
   # psych::KMO(resid_corb)$MSA
   # psych::fa.parallel(resid_corb, fm = "pa", fa = "fa",
   #                    n.obs = lavTech(msat_fit, "nclusters")[[1]], 
   #                    n.iter = 30 * nitem, 
   #                    plot = FALSE)$nfact
   # Dimensionality
   corw <- lavTech(msat_fit, what = "cor.ov")$within
   corb <- lavTech(msat_fit, what = "cor.ov")$id
   paw <- do.call(fa.parallel, 
                  args = c(list(x = corw, 
                                n.obs = nobs(msat_fit) -
                                   lavTech(msat_fit, "nclusters")[[1]],
                                n.iter = 30 * nitem,
                                plot = FALSE), pa_args))
   pab <- do.call(fa.parallel, 
                  args = c(list(x = corw, 
                                n.obs = lavTech(msat_fit, "nclusters")[[1]],
                                n.iter = 30 * nitem,
                                plot = FALSE), pa_args))
   if (pa_args$fa == "pc") {
      ncompw <- paw$ncomp
      ncompb <- pab$ncomp
   } else if (pa_args$fa == "fa") {
      ncompw <- paw$nfact
      ncompb <- pab$nfact
   }
   list(alpha = c(alpha2l = alpha2l, alphab = alphab, alphaw = alphaw), 
        alpha_ci = do.call(rbind, sim_alpha_cis), 
        omega = c(omega2l = omega2l, omegab = omegab, omegaw = omegaw), 
        omega_ci = do.call(rbind, sim_omega_cis), 
        ncomp = c(within = ncompw, between = ncompb))
}

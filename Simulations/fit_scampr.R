## Function fit LGCP via mgcv

fit_scampr <- function(data, quad) {
  
  # set up the PO and quadrature as a single data frame
  quad$present <- 0
  tmp.dat <- rbind(data[,c("x", "y", "env", "bias", "quad.size", "present")], quad[,c("x", "y", "env", "bias", "quad.size", "present")])

  # fit the base model without SRE
  m0 <- scampr(formula = present ~ env, data = tmp.dat, include.sre = F)
  # perform the basis opt.
  m_va <- basis.search.po(m0, domain.data = quad, return.model = T)
  m_lp <- basis.search.po(m0, domain.data = quad, return.model = T, which.approx = "laplace")
  
  ## EVALUATION OF MODEL #######################################################
  
  # get the predictions
  m.prd_va <- exp(m_va$fitted.values[tmp.dat$present == 0])
  m.prd_lp <- exp(m_lp$fitted.values[tmp.dat$present == 0])
  # calculate the KL divergence metric
  KLdiv_va <- as.numeric(quad$quad.size %*% (quad$lambda * log(quad$lambda / m.prd_va))) - as.numeric(quad$quad.size %*% (quad$lambda - m.prd_va))
  KLdiv_lp <- as.numeric(quad$quad.size %*% (quad$lambda * log(quad$lambda / m.prd_lp))) - as.numeric(quad$quad.size %*% (quad$lambda - m.prd_lp))
  # calculate the relative MAE
  MAE_va <- mean(abs(m.prd_va-quad$lambda))
  MAE_lp <- mean(abs(m.prd_lp-quad$lambda))
  # get the fixed par estimate
  BETA_HAT_va <- unname(m_va$coefficients["env"])
  BETA_HAT_lp <- unname(m_lp$coefficients["env"])
  # and the estimate's SE
  BETA_SE_va <- m_va$fixed.effects["env", "Std. Error"]
  BETA_SE_lp <- m_lp$fixed.effects["env", "Std. Error"]
  # calculate the rmse in beta estimate
  SQER_BETA_va <- (BETA_HAT_va - attr(pres, "sim.info")$Beta.env)^2
  SQER_BETA_lp <- (BETA_HAT_lp - attr(pres, "sim.info")$Beta.env)^2
  # calculate coverage for beta estimate
  COVER_BETA_va <- BETA_HAT_va + qnorm(0.025) * BETA_SE_va  <= attr(pres, "sim.info")$Beta.env & BETA_HAT_va + qnorm(0.975) * BETA_SE_va >= attr(pres, "sim.info")$Beta.env
  COVER_BETA_lp <- BETA_HAT_lp + qnorm(0.025) * BETA_SE_lp  <= attr(pres, "sim.info")$Beta.env & BETA_HAT_lp + qnorm(0.975) * BETA_SE_lp >= attr(pres, "sim.info")$Beta.env
  # calculate the rmse in range parameter estimate
  SQER_RHO <- c(NA, NA)
  # calculate coverage for range parameter estimate
  COVER_RHO <- c(NA, NA)
  e0 <- data.frame(GP_APPROX = c("OPT", "OPT"), COV_FN = c(NA, NA), POW = c(NA, NA),
                   K = c(NA, NA), FIT = rep("scampr", 2), KL = c(KLdiv_va, KLdiv_lp), MAE = c(MAE_va, MAE_lp),
                   SQER_BETA = c(SQER_BETA_va, SQER_BETA_lp), COVER_BETA = c(COVER_BETA_va, COVER_BETA_lp), BETA_HAT = c(BETA_HAT_va, BETA_HAT_lp),
                   SQER_RHO = c(NA, NA), COVER_RHO = c(NA, NA), RHO_HAT = c(NA, NA),
                   TIME = c(m_va$cpu["basis.search"], m_lp$cpu["basis.search"]), DEF_KNOT_LOC = c(T, T),
                   INLA_K = c(if (is.null(m_va$basis.functions)) {0} else {nrow(m_va$basis.functions)}, if (is.null(m_lp$basis.functions)) {0} else {nrow(m_lp$basis.functions)}),
                   CRIT = c(logLik(m_va), logLik(m_lp)), EDF = c(NA, NA), method = c("VA", "LP")
  )

  return(e0)
}
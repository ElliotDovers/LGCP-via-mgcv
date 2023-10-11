## Function fit LGCP via mgcv

fit_gam <- function(data, quad, pred, n_knots, which.cov.fn = 2, crit.method = c("REML", "GCV.Cp")) {
  
  crit.method <- match.arg(crit.method)
  
  # set up the PO and quadrature as a single data frame
  quad$present <- 0
  tmp.dat <- rbind(data[,c("x", "y", "env", "bias", "quad.size", "present")], quad[,c("x", "y", "env", "bias", "quad.size", "present")])
  # calculate the point weights for DWPR
  tmp.dat$p.wt <- tmp.dat$quad.size
  tmp.dat$p.wt[tmp.dat$present == 1] <- 1e-6
  # set up the adjusted response
  tmp.dat$response <- tmp.dat$present / tmp.dat$p.wt # retaining this approach as the speed is very slightly better
  # check the IPP fit
  # mA <- gam(response ~ env, data = tmp.dat, family = poisson(), weights = p.wt, method="REML")
  # mB <- scampr::scampr(present ~ env, data = tmp.dat, include.sre = F)
  time0 <- system.time(assign("m0", gam(response ~ env + te(x,y, bs = "gp", k = n_knots, m = list(c(which.cov.fn)), d = 2), data = tmp.dat, family = poisson(), weights = p.wt, method=crit.method)))

  ## EVALUATION OF MODEL #######################################################
  
  # get the predictions
  m.prd <- m0$fitted.values[tmp.dat$present == 0]
  # calculate the KL divergence metric
  KLdiv <- as.numeric(pred$quad.size %*% (pred$lambda * log(pred$lambda / m.prd))) - as.numeric(pred$quad.size %*% (pred$lambda - m.prd))
  # calculate the relative MAE
  MAE <- mean(abs(m.prd-pred$lambda))
  # get the fixed par estimate
  BETA_HAT <- unname(m0$coefficients["env"])
  # and the estimate's SE
  BETA_SE <- sqrt(diag(vcov(m0, unconditional = T)))["env"]
  # calculate the rmse in beta estimate
  SQER_BETA <- (BETA_HAT - attr(pres, "sim.info")$Beta.env)^2
  # calculate coverage for beta estimate
  COVER_BETA <- BETA_HAT + qnorm(0.025) * BETA_SE  <= attr(pres, "sim.info")$Beta.env & BETA_HAT + qnorm(0.975) * BETA_SE >= attr(pres, "sim.info")$Beta.env
  # calculate the rmse in range parameter estimate
  SQER_RHO <- (m0$smooth[[1]]$margin[[1]]$gp.defn[2] - attr(pres, "sim.info")$latent.practical.range)^2
  # calculate coverage for range parameter estimate
  COVER_RHO <- NA
  e0 <- data.frame(GP_APPROX = class(m0$smooth[[1]]$margin[[1]]), COV_FN = m0$smooth[[1]]$margin[[1]]$gp.defn[1], POW = m0$smooth[[1]]$margin[[1]]$gp.defn[3],
                   K = m0$smooth[[1]]$margin[[1]]$bs.dim, FIT = "mgcv", KL = KLdiv, MAE = MAE,
                   SQER_BETA = SQER_BETA, COVER_BETA = COVER_BETA, BETA_HAT = BETA_HAT,
                   SQER_RHO = SQER_RHO, COVER_RHO = COVER_RHO, RHO_HAT = m0$smooth[[1]]$margin[[1]]$gp.defn[2],
                   TIME = time0[3], DEF_KNOT_LOC = T, INLA_K = NA, CRIT = NA, EDF = sum(m0$edf), method = crit.method
  )

  return(e0)
}
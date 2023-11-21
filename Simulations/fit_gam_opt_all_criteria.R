## Function fit LGCP via mgcv

fit_gam_opt <- function(data, quad, pred, n_knots, which.cov.fn = 2, crit.method = c("REML", "GCV.Cp"), tolerance = 2){
  
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
  # if (missing(n_knots)) {
  #   n_knots <- 25
  # }
  # n_knots <<- n_knots
  
  # set the functions to be optimised
  fit.new.ll <- function(logrho) {
    tmp.m <- gam(response ~ env + te(x,y, bs = "gp", k = n_knots, m = list(c(which.cov.fn, exp(logrho), 1)), d = 2), data=tmp.dat, family=poisson(), weights=p.wt, method=crit.method)
    return(-logLik(tmp.m))
  }
  fit.new.aic <- function(logrho) {
    tmp.m <- gam(response ~ env + te(x,y, bs = "gp", k = n_knots, m = list(c(which.cov.fn, exp(logrho), 1)), d = 2), data=tmp.dat, family=poisson(), weights=p.wt, method=crit.method)
    return(AIC(tmp.m))
  }
  fit.new.score <- function(logrho) {
    tmp.m <- gam(response ~ env + te(x,y, bs = "gp", k = n_knots, m = list(c(which.cov.fn, exp(logrho), 1)), d = 2), data=tmp.dat, family=poisson(), weights=p.wt, method=crit.method)
    return(tmp.m$gcv.ubre)
  }
  fit.new.bic <- function(logrho) {
    tmp.m <- gam(response ~ env + te(x,y, bs = "gp", k = n_knots, m = list(c(which.cov.fn, exp(logrho), 1)), d = 2), data=tmp.dat, family=poisson(), weights=p.wt, method=crit.method)
    return(AIC(tmp.m, k = log(sum(tmp.dat$present)))) # adjusted to be proper BIC (using number of presence points)
  }
  # optimise the range parameter
  time0 <- system.time(assign("opt0", optimize(fit.new.ll, interval = log(c(0.1, 140)), tol = tolerance)))
  time1 <- system.time(assign("opt1", optimize(fit.new.aic, interval = log(c(0.1, 140)), tol = tolerance)))
  time2 <- system.time(assign("opt2", optimize(fit.new.score, interval = log(c(0.1, 140)), tol = tolerance)))
  time3 <- system.time(assign("opt3", optimize(fit.new.bic, interval = log(c(0.1, 140)), tol = tolerance)))
  # fit the models
  time0B <- system.time(assign("m0", gam(response ~ env + te(x,y, bs = "gp", k = n_knots, m = list(c(which.cov.fn, exp(opt0$minimum), 1)), d = 2), data = tmp.dat, family = poisson(), weights = p.wt, method=crit.method)))
  time1B <- system.time(assign("m1", gam(response ~ env + te(x,y, bs = "gp", k = n_knots, m = list(c(which.cov.fn, exp(opt1$minimum), 1)), d = 2), data = tmp.dat, family = poisson(), weights = p.wt, method=crit.method)))
  time2B <- system.time(assign("m2", gam(response ~ env + te(x,y, bs = "gp", k = n_knots, m = list(c(which.cov.fn, exp(opt2$minimum), 1)), d = 2), data = tmp.dat, family = poisson(), weights = p.wt, method=crit.method)))
  time3B <- system.time(assign("m3", gam(response ~ env + te(x,y, bs = "gp", k = n_knots, m = list(c(which.cov.fn, exp(opt3$minimum), 1)), d = 2), data = tmp.dat, family = poisson(), weights = p.wt, method=crit.method)))
  
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
                   TIME = time0[3] + time0B[3], DEF_KNOT_LOC = T, INLA_K = NA, CRIT = "ll", EDF = sum(m0$edf), method = crit.method
  )
  
  # get the predictions
  m.prd <- m1$fitted.values[tmp.dat$present == 0]
  # calculate the KL divergence metric
  KLdiv <- as.numeric(pred$quad.size %*% (pred$lambda * log(pred$lambda / m.prd))) - as.numeric(pred$quad.size %*% (pred$lambda - m.prd))
  # calculate the relative MAE
  MAE <- mean(abs(m.prd-pred$lambda))
  # get the fixed par estimate
  BETA_HAT <- unname(m1$coefficients["env"])
  # and the estimate's SE
  BETA_SE <- sqrt(diag(vcov(m1, unconditional = T)))["env"]
  # calculate the rmse in beta estimate
  SQER_BETA <- (BETA_HAT - attr(pres, "sim.info")$Beta.env)^2
  # calculate coverage for beta estimate
  COVER_BETA <- BETA_HAT + qnorm(0.025) * BETA_SE  <= attr(pres, "sim.info")$Beta.env & BETA_HAT + qnorm(0.975) * BETA_SE >= attr(pres, "sim.info")$Beta.env
  # calculate the rmse in range parameter estimate
  SQER_RHO <- (m1$smooth[[1]]$margin[[1]]$gp.defn[2] - attr(pres, "sim.info")$latent.practical.range)^2
  # calculate coverage for range parameter estimate
  COVER_RHO <- NA
  e1 <- data.frame(GP_APPROX = class(m1$smooth[[1]]$margin[[1]]), COV_FN = m1$smooth[[1]]$margin[[1]]$gp.defn[1], POW = m1$smooth[[1]]$margin[[1]]$gp.defn[3],
                   K = m1$smooth[[1]]$margin[[1]]$bs.dim, FIT = "mgcv", KL = KLdiv, MAE = MAE,
                   SQER_BETA = SQER_BETA, COVER_BETA = COVER_BETA, BETA_HAT = BETA_HAT,
                   SQER_RHO = SQER_RHO, COVER_RHO = COVER_RHO, RHO_HAT = m1$smooth[[1]]$margin[[1]]$gp.defn[2],
                   TIME = time1[3] + time1B[3], DEF_KNOT_LOC = T, INLA_K = NA, CRIT = "aic", EDF = sum(m1$edf), method = crit.method
  )
  
  # get the predictions
  m.prd <- m2$fitted.values[tmp.dat$present == 0]
  # calculate the KL divergence metric
  KLdiv <- as.numeric(pred$quad.size %*% (pred$lambda * log(pred$lambda / m.prd))) - as.numeric(pred$quad.size %*% (pred$lambda - m.prd))
  # calculate the relative MAE
  MAE <- mean(abs(m.prd-pred$lambda))
  # get the fixed par estimate
  BETA_HAT <- unname(m2$coefficients["env"])
  # and the estimate's SE
  BETA_SE <- sqrt(diag(vcov(m2, unconditional = T)))["env"]
  # calculate the rmse in beta estimate
  SQER_BETA <- (BETA_HAT - attr(pres, "sim.info")$Beta.env)^2
  # calculate coverage for beta estimate
  COVER_BETA <- BETA_HAT + qnorm(0.025) * BETA_SE  <= attr(pres, "sim.info")$Beta.env & BETA_HAT + qnorm(0.975) * BETA_SE >= attr(pres, "sim.info")$Beta.env
  # calculate the rmse in range parameter estimate
  SQER_RHO <- (m2$smooth[[1]]$margin[[1]]$gp.defn[2] - attr(pres, "sim.info")$latent.practical.range)^2
  # calculate coverage for range parameter estimate
  COVER_RHO <- NA
  e2 <- data.frame(GP_APPROX = class(m2$smooth[[1]]$margin[[1]]), COV_FN = m2$smooth[[1]]$margin[[1]]$gp.defn[1], POW = m2$smooth[[1]]$margin[[1]]$gp.defn[3],
                   K = m2$smooth[[1]]$margin[[1]]$bs.dim, FIT = "mgcv", KL = KLdiv, MAE = MAE,
                   SQER_BETA = SQER_BETA, COVER_BETA = COVER_BETA, BETA_HAT = BETA_HAT,
                   SQER_RHO = SQER_RHO, COVER_RHO = COVER_RHO, RHO_HAT = m2$smooth[[1]]$margin[[1]]$gp.defn[2],
                   TIME = time2[3] + time2B[3], DEF_KNOT_LOC = T, INLA_K = NA, CRIT = "score", EDF = sum(m2$edf), method = crit.method
  )
  
  # get the predictions
  m.prd <- m3$fitted.values[tmp.dat$present == 0]
  # calculate the KL divergence metric
  KLdiv <- as.numeric(pred$quad.size %*% (pred$lambda * log(pred$lambda / m.prd))) - as.numeric(pred$quad.size %*% (pred$lambda - m.prd))
  # calculate the relative MAE
  MAE <- mean(abs(m.prd-pred$lambda))
  # get the fixed par estimate
  BETA_HAT <- unname(m3$coefficients["env"])
  # and the estimate's SE
  BETA_SE <- sqrt(diag(vcov(m3, unconditional = T)))["env"]
  # calculate the rmse in beta estimate
  SQER_BETA <- (BETA_HAT - attr(pres, "sim.info")$Beta.env)^2
  # calculate coverage for beta estimate
  COVER_BETA <- BETA_HAT + qnorm(0.025) * BETA_SE  <= attr(pres, "sim.info")$Beta.env & BETA_HAT + qnorm(0.975) * BETA_SE >= attr(pres, "sim.info")$Beta.env
  # calculate the rmse in range parameter estimate
  SQER_RHO <- (m3$smooth[[1]]$margin[[1]]$gp.defn[2] - attr(pres, "sim.info")$latent.practical.range)^2
  # calculate coverage for range parameter estimate
  COVER_RHO <- NA
  e3 <- data.frame(GP_APPROX = class(m3$smooth[[1]]$margin[[1]]), COV_FN = m3$smooth[[1]]$margin[[1]]$gp.defn[1], POW = m3$smooth[[1]]$margin[[1]]$gp.defn[3],
                   K = m3$smooth[[1]]$margin[[1]]$bs.dim, FIT = "mgcv", KL = KLdiv, MAE = MAE,
                   SQER_BETA = SQER_BETA, COVER_BETA = COVER_BETA, BETA_HAT = BETA_HAT,
                   SQER_RHO = SQER_RHO, COVER_RHO = COVER_RHO, RHO_HAT = m3$smooth[[1]]$margin[[1]]$gp.defn[2],
                   TIME = time2[3] + time2B[3], DEF_KNOT_LOC = T, INLA_K = NA, CRIT = "bic", EDF = sum(m3$edf), method = crit.method
  )

  ret.obj <- rbind(e0, e1, e2, e3)
  return(ret.obj)
}
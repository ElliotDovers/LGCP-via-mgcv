## Function to fit LGCP via INLA (two ways)

fit_inla <- function(data, quad, n_mesh, pred){
  
  # set the mesh according to n_mesh
  mesh <- inla.mesh.2d(loc.domain = quad[ , c("x", "y")], max.edge=c(5,10), cutoff=2, offset = c(5,10), max.n.strict = c(n_mesh, 0))
  
  # if (missing(n_mesh)) {
  #   mesh <- inla.mesh.2d(loc.domain = quad[ , c("x", "y")], max.edge=c(10,30), cutoff=2, offset = c(5,20))
  #   n_mesh <- mesh$n
  #   def_knot_loc <- T
  # } else {
  #   knots <- data.frame(x = runif(n_mesh, 0, 100), y = runif(n_mesh, 0, 100))
  #   mesh <- inla.mesh.create(loc = knots)
  #   def_knot_loc <- F
  # }
  
  # set the spde representation to be the mesh with penalised complexity priors
  spde.pcmatern <- inla.spde2.pcmatern(mesh,
                                       alpha = 1.5,
                                       prior.sigma = c(10, 0.01),
                                       prior.range = c(1, 0.01)
  )
  # also try the default settings
  spde.matern <- inla.spde2.matern(mesh)
  
  # make A matrix for unstructured data
  data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(data[,c("x","y")]))
  
  # make A matrix for quadrature points (NOTE THIS WILL BE DIAGONAL IF MESH == QUAD)
  quad_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(quad[,c("x","y")]))
  
  # make A matrix for the prediction points
  pred_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(pred[ , c("x","y")]))
  
  # LGCP
  
  # One spatial field
  # Uses Simpson approach for PP data
  
  # set the number of integration points, presence points in PO data and prediction points
  nq <- nrow(quad)
  n <- nrow(data)
  np <- nrow(pred)
  
  # change data to include 0s for nodes and 1s for presences
  y.pp <- rep(0:1, c(nq, n))
  
  # add expectation vector (area for integration points/nodes and 0 for presences)
  e.pp <- c(quad$quad.size, rep(0, n))
  
  # combine integration point A matrix over quadrature with PO data A matrix
  A.pp <- rbind(quad_data_A, data_A)
  
  # Create data stack
  stk_data <- inla.stack(data=list(y=y.pp, e = e.pp),
                         effects=list(list(data.frame(Intercept=rep(1,nq+n)), env = c(quad$env, data$env)), list(i=1:mesh$n)),
                         A=list(1,A.pp),
                         tag="po_data")
  
  # create the prediction stack
  stk_pred_response <- inla.stack(data=list(y=NA),
                                  effects = list(list(data.frame(Intercept=rep(1,np))), env = pred$env, list(i=1:mesh$n)),
                                  A=list(1,1, pred_data_A),
                                  tag='pred_response')
  
  # combine the stacks
  stk <- inla.stack(stk_data, stk_pred_response)
  
  # fit the model
  result <- inla(y ~ Intercept + env + f(i, model = spde.pcmatern) - 1,
                 family="poisson",
                 data=inla.stack.data(stk),
                 control.predictor=list(A=inla.stack.A(stk), compute=TRUE),
                 control.family = list(link = "log"),
                 E = inla.stack.data(stk)$e,
                 control.compute = list(cpo=TRUE, waic = TRUE, dic = TRUE)
  )
  resultB <- inla(y ~ Intercept + env + f(i, model = spde.matern) - 1,
                  family="poisson",
                  data=inla.stack.data(stk),
                  control.predictor=list(A=inla.stack.A(stk), compute=TRUE),
                  control.family = list(link = "log"),
                  E = inla.stack.data(stk)$e,
                  control.compute = list(cpo=TRUE, waic = TRUE, dic = TRUE)
  )
  
  # create index to extract predictions
  index.pred.response <- inla.stack.index(stk, tag="pred_response")$data
  # get the predictions
  m.prd <- exp(result$summary.fitted.values$mean[index.pred.response])
  # calculate the KL divergence metric
  KLdiv <- as.numeric(pred$quad.size %*% (pred$lambda * log(pred$lambda / m.prd))) - as.numeric(pred$quad.size %*% (pred$lambda - m.prd))
  # calculate the relative MAE
  MAE <- mean(abs(m.prd-pred$lambda))
  # get the fixed par estimate
  BETA_HAT <- result$summary.fixed["env","mean"]
  # and the estimate's SE
  BETA_SE <- result$summary.fixed["env","sd"] # here we get the standard dev. of distribution of Beta so can use that
  # calculate the rmse in beta estimate
  SQER_BETA <- (BETA_HAT - attr(pres, "sim.info")$Beta.env)^2
  # calculate coverage for beta estimate (credible)
  # COVER_BETA <- result$summary.fixed["env","0.025quant"] <= attr(pres, "sim.info")$Beta.env & result$summary.fixed["env","0.975quant"] >= attr(pres, "sim.info")$Beta.env
  # calculate coverage for beta estimate (WALD)
  COVER_BETA <- BETA_HAT + qnorm(0.025) * BETA_SE  <= attr(pres, "sim.info")$Beta.env & BETA_HAT + qnorm(0.975) * BETA_SE >= attr(pres, "sim.info")$Beta.env
  # calculate the rmse in range parameter estimate
  SQER_RHO <- (result$summary.hyperpar["Range","mean"] - attr(pres, "sim.info")$latent.practical.range)^2
  # calculate coverage for range parameter estimate
  COVER_RHO <- result$summary.hyperpar["Range","0.025quant"]/2 <= attr(pres, "sim.info")$latent.practical.range & result$summary.hyperpar["Range","0.975quant"]/2 >= attr(pres, "sim.info")$latent.practical.range
  
  e0 <- data.frame(GP_APPROX = "pcmatern", COV_FN = 2, POW = 1, K = n_mesh, FIT = "inla", KL = KLdiv, MAE = MAE,
                   SQER_BETA = SQER_BETA, COVER_BETA = COVER_BETA, BETA_HAT = BETA_HAT,
                   SQER_RHO = SQER_RHO, COVER_RHO = COVER_RHO, RHO_HAT = result$summary.hyperpar["Range","mean"],
                   TIME = result$cpu.used[4], DEF_KNOT_LOC = T, INLA_K = mesh$n, CRIT = NA, EDF = NA, method = NA
  )
  
  # get the predictions
  m.prd <- exp(resultB$summary.fitted.values$mean[index.pred.response])
  # calculate the KL divergence metric
  KLdiv <- as.numeric(pred$quad.size %*% (pred$lambda * log(pred$lambda / m.prd))) - as.numeric(pred$quad.size %*% (pred$lambda - m.prd))
  # calculate the relative MAE
  MAE <- mean(abs(m.prd-pred$lambda))
  # get the fixed par estimate
  BETA_HAT <- resultB$summary.fixed["env","mean"]
  # calculate the rmse in beta estimate
  SQER_BETA <- (BETA_HAT - attr(pres, "sim.info")$Beta.env)^2
  # calculate coverage for beta estimate
  COVER_BETA <- resultB$summary.fixed["env","0.025quant"] <= attr(pres, "sim.info")$Beta.env & resultB$summary.fixed["env","0.975quant"] >= attr(pres, "sim.info")$Beta.env
  # calculate the rmse in range parameter estimate
  SQER_RHO <- (sqrt(8)/exp(resultB$summary.hyperpar["Theta2","mean"]) - attr(pres, "sim.info")$latent.practical.range)^2
  # calculate coverage for range parameter estimate (note the flip of the conditions due to the inversion!)
  COVER_RHO <- (sqrt(8)/exp(resultB$summary.hyperpar["Theta2","0.025quant"]))/2 >= attr(pres, "sim.info")$latent.practical.range & (sqrt(8)/exp(resultB$summary.hyperpar["Theta2","0.975quant"]))/2 <= attr(pres, "sim.info")$latent.practical.range
  
  e1 <- data.frame(GP_APPROX = "default", COV_FN = 2, POW = 1, K = n_mesh, FIT = "inla", KL = KLdiv, MAE = MAE,
                   SQER_BETA = SQER_BETA, COVER_BETA = COVER_BETA, BETA_HAT = BETA_HAT,
                   SQER_RHO = SQER_RHO, COVER_RHO = COVER_RHO, RHO_HAT = (sqrt(8)/exp(resultB$summary.hyperpar["Theta2","mean"])),
                   TIME = resultB$cpu.used[4], DEF_KNOT_LOC = T, INLA_K = mesh$n, CRIT = NA, EDF = NA, method = NA
  )
  
  return(rbind(e0, e1))
}
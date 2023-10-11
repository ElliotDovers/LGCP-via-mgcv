library(mgcv)
library(RandomFields)
library(spatstat)
library(sp)
library(INLA)
library(fields)

################################################################################
# Function to interpolate some covariate at x, y locations #####################
################################################################################
interp.covar <- function(x.loc, y.loc, covar.name, domain.data){
  
  # turn the quadrature into a spatial pixels data frame
  sp.domain <- sp::SpatialPixelsDataFrame(points = domain.data[,c("x", "y")], data = domain.data[ , !colnames(domain.data) %in% c("x", "y", "quad.size")])
  
  # turn coordinates into SpatialPoints object:
  spp = sp::SpatialPoints(data.frame(x = x.loc,y = y.loc)) 
  
  # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
  v <- sp::over(spp, sp.domain[ , covar.name])
  v[is.na(v)] = 0 # NAs are a problem! Remove them
  return(v[,1])
}
################################################################################

# Get the job array
tab <- read.csv("job_array_all.csv")

## THE JOB INDEX HAS BEEN FIXED FOR DEMONSTRATIVE PURPOSES
job = 1
# # determine job number from pbs script
# job = as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

################################################################################
# Set parameters that define the scenarios:

# random seed / simulation number
seed = tab$sim[tab$job == job]
# environment range of effect
env.range = tab$env_range[tab$job == job]
# latent range of effect
lat.range = tab$lat_range[tab$job == job]
# model
model_to_test <- tab$fit_model[tab$job == job]
# model
beta0 <- tab$intercept[tab$job == job]

source("sim_lgcp_pp.R")
pres <- sim_lgcp_pp(Intercept = beta0,
                    rseed = seed,
                    env.covariate.type = "random_field",
                    env.covariate.range = env.range,
                    env.covariate.covar.function = "stable",
                    latent.field = T,
                    latent.practical.range = lat.range,
                    latent.smoothness = 0.5, # fixed to be an exponential covariance function
                    plotting = F
)
# set the entire domain truth grid to be the quadrature and the prediction points
quad <- pred <- attr(pres, "truth.grid")

# # set up the inla meshes
# mesh0 <- inla.mesh.2d(loc.domain = quad[ , c("x", "y")], max.edge=c(10,30), cutoff=2, offset = c(5,20))
# mesh1 <- inla.mesh.create(loc = data.frame(x = runif(100, 0, 100), y = runif(100, 0, 100)))
# mesh2 <- inla.mesh.create(loc = data.frame(x = runif(200, 0, 100), y = runif(200, 0, 100)))
# mesh3 <- inla.mesh.create(loc = data.frame(x = runif(300, 0, 100), y = runif(300, 0, 100)))
# mesh4 <- inla.mesh.create(loc = data.frame(x = runif(400, 0, 100), y = runif(400, 0, 100)))
# sum(in.out(matrix(c(0,0,100,0,100,100,0,100,0,0), ncol = 2, byrow = T), mesh0$loc[,1:2])) # 260

# fit and predict using the appropriate model for the job
if (model_to_test == "inla") {
  source("fit_inla.R")
  res0 <- fit_inla(data = pres, quad = quad, n_mesh = 25, pred = pred)
  res1 <- fit_inla(data = pres, quad = quad, n_mesh = 100, pred = pred)
  res2 <- fit_inla(data = pres, quad = quad, n_mesh = 200, pred = pred)
  res3 <- fit_inla(data = pres, quad = quad, n_mesh = 300, pred = pred)
  res4 <- fit_inla(data = pres, quad = quad, n_mesh = 400, pred = pred)
  # combine results
  res <- rbind(res0, res1, res2, res3, res4)
} else if (model_to_test == "gam_opt") {
  source("fit_gam_opt.R")
  res0 <- fit_gam_opt(data = pres, quad = quad, n_knots = 25, pred = pred)
  res1 <- fit_gam_opt(data = pres, quad = quad, n_knots = 100, pred = pred)
  res2 <- fit_gam_opt(data = pres, quad = quad, n_knots = 200, pred = pred)
  res3 <- fit_gam_opt(data = pres, quad = quad, n_knots = 300, pred = pred)
  res4 <- fit_gam_opt(data = pres, quad = quad, n_knots = 400, pred = pred)
  res5 <- fit_gam_opt(data = pres, quad = quad, n_knots = 25, pred = pred, crit.method = "GCV.Cp")
  res6 <- fit_gam_opt(data = pres, quad = quad, n_knots = 100, pred = pred, crit.method = "GCV.Cp")
  res7 <- fit_gam_opt(data = pres, quad = quad, n_knots = 200, pred = pred, crit.method = "GCV.Cp")
  res8 <- fit_gam_opt(data = pres, quad = quad, n_knots = 300, pred = pred, crit.method = "GCV.Cp")
  res9 <- fit_gam_opt(data = pres, quad = quad, n_knots = 400, pred = pred, crit.method = "GCV.Cp")
  # combine results
  res <- rbind(res0, res1, res2, res3, res4, res5, res6, res7, res8, res9)
} else if (model_to_test == "mgcv") {
  source("fit_gam.R")
  res0 <- fit_gam(data = pres, quad = quad, n_knots = 25, pred = pred)
  res1 <- fit_gam(data = pres, quad = quad, n_knots = 100, pred = pred)
  res2 <- fit_gam(data = pres, quad = quad, n_knots = 200, pred = pred)
  res3 <- fit_gam(data = pres, quad = quad, n_knots = 300, pred = pred)
  res4 <- fit_gam(data = pres, quad = quad, n_knots = 400, pred = pred)
  res5 <- fit_gam(data = pres, quad = quad, n_knots = 25, pred = pred, crit.method = "GCV.Cp")
  res6 <- fit_gam(data = pres, quad = quad, n_knots = 100, pred = pred, crit.method = "GCV.Cp")
  res7 <- fit_gam(data = pres, quad = quad, n_knots = 200, pred = pred, crit.method = "GCV.Cp")
  res8 <- fit_gam(data = pres, quad = quad, n_knots = 300, pred = pred, crit.method = "GCV.Cp")
  res9 <- fit_gam(data = pres, quad = quad, n_knots = 400, pred = pred, crit.method = "GCV.Cp")
  # combine results
  res <- rbind(res0, res1, res2, res3, res4, res5, res6, res7, res8, res9)
} else if (model_to_test == "youngman") {
  source("fit_youngman.R")
  res0 <- fit_youngman(data = pres, quad = quad, n_knots = 25, pred = pred)
  res1 <- fit_youngman(data = pres, quad = quad, n_knots = 100, pred = pred)
  res2 <- fit_youngman(data = pres, quad = quad, n_knots = 200, pred = pred)
  res3 <- fit_youngman(data = pres, quad = quad, n_knots = 300, pred = pred)
  res4 <- fit_youngman(data = pres, quad = quad, n_knots = 400, pred = pred)
  res5 <- fit_youngman(data = pres, quad = quad, n_knots = 25, pred = pred, crit.method = "GCV.Cp")
  res6 <- fit_youngman(data = pres, quad = quad, n_knots = 100, pred = pred, crit.method = "GCV.Cp")
  res7 <- fit_youngman(data = pres, quad = quad, n_knots = 200, pred = pred, crit.method = "GCV.Cp")
  res8 <- fit_youngman(data = pres, quad = quad, n_knots = 300, pred = pred, crit.method = "GCV.Cp")
  res9 <- fit_youngman(data = pres, quad = quad, n_knots = 400, pred = pred, crit.method = "GCV.Cp")
  # combine results
  res <- rbind(res0, res1, res2, res3, res4, res5, res6, res7, res8, res9)
}

# collate with sim info
res_tab <- cbind(tab[tab$job == job, ], res)
res_tab$n <- nrow(pres)

# save the simulation result table in folder "Results" inside the base dir
save(list = "res_tab", file = paste0(getwd(), "/Results/res_", job, ".RDATA"))

# set the root working directory as the current one
home.wd <- getwd()
plot.res <- 500

# install the required software packages #######################################
if(!require(mgcv, quietly = T)){
  install.packages("mgcv")
  library(mgcv)
}
if(!require(spatstat, quietly = T)){
  install.packages("spatstat")
  library(spatstat)
}

data(gorillas, package = "spatstat.data")
pp = data.frame(gorillas,
                lapply(gorillas.extra,function(x){x[gorillas]}),pt=1,wt=1e-6)
q_xy = data.frame(gorillas.extra[[1]])[,c("x","y")] # extract x and y
quad = data.frame(q_xy,lapply(gorillas.extra,function(x){x[q_xy]}),
                  pt=0,wt=area(gorillas$window)/nrow(q_xy))
dat = merge(pp, quad, all=T)

# center and scale covariates
dat$elevation <- scale(dat$elevation)
dat$slopeangle <- scale(dat$slopeangle)
dat$waterdist <- scale(dat$waterdist)

# fit the models with two candidate basis dimensions
# m_k200_te <- gam(pt/wt ~ elevation + waterdist + slopeangle + heat + slopetype + vegetation + te(x, y, bs = "gp", k = 200, m = list(3), d = 2), data=dat, family=poisson(), weights=wt, method="REML")
# m_k400_te <- gam(pt/wt ~ elevation + waterdist + slopeangle + heat + slopetype + vegetation + te(x, y, bs = "gp", k = 400, m = list(3), d = 2), data=dat, family=poisson(), weights=wt, method="REML")
m_k200 <- gam(pt/wt ~ elevation + waterdist + slopeangle + heat + slopetype + vegetation + s(x, y, bs = "gp", k = 200, m = 3), data=dat, family=poisson(), weights=wt, method="REML")
m_k400 <- gam(pt/wt ~ elevation + waterdist + slopeangle + heat + slopetype + vegetation + s(x, y, bs = "gp", k = 400, m = 3), data=dat, family=poisson(), weights=wt, method="REML")

# compare the edf to k ratios
sum(m_k200$edf) / 200
sum(m_k400$edf) / 400

dists <- nndist(gorillas)
range_interval <- range(dists[dists != 0])
# set up the function to be minimized
objective_fn = function(rho) {
  tmp.m = gam(pt/wt ~ elevation + waterdist + slopeangle + heat + slopetype + vegetation +
                s(x, y, bs = "gp", k = 400, m = c(3,rho)),
              data=dat, family=poisson(), weights=wt, method="REML"
  )
  return(tmp.m$gcv.ubre) # the "method" specific criterion
}
# find the optimized range parameter
if (file.exists("pre-calculated_optim.RDATA")) {
  load("pre-calculated_optim.RDATA")
} else {
  opt_time <- system.time(assign("opt", optimize(objective_fn, interval = range_interval)))
  save(list = c("opt", "opt_time"), file = "pre-calculated_optim.RDATA")
}

# m <- gam(pt/wt ~ elevation + waterdist + slopeangle + heat + slopetype + vegetation + te(x, y, bs = "gp", k = 400, m = list(c(3, opt$minimum)), d = 2), data=dat, family=poisson(), weights=wt, method="REML")
m <- gam(pt/wt ~ elevation + waterdist + slopeangle + heat + slopetype + vegetation + s(x, y, bs = "gp", k = 400, m = c(3, opt$minimum)), data=dat, family=poisson(), weights=wt, method="REML")

# fit an IPP model to contrast with the fitted LGCP
m_ipp <- gam(pt/wt ~ elevation + waterdist + slopeangle + heat + slopetype + vegetation , data=dat, family=poisson(), weights=wt, method="REML")

# set the domain data points (in this case the quadrature we used)
domain.grid <- dat[dat$pt == 0, ]

# predict intensity values
domain.grid$z_ipp = predict(m_ipp, newdata=domain.grid, type = "response")
domain.grid$z = predict(m, newdata=domain.grid, type = "response")

# create the pixel images (uses the window supplied in the original gorillas data)
pred_ipp.im = as.im(domain.grid[,c("x","y","z_ipp")], W = gorillas$window)
pred.im = as.im(domain.grid[,c("x","y","z")], W = gorillas$window)

# calculate the observed K functions
K_obs_ipp <- Kinhom(gorillas, lambda = pred_ipp.im, correction = "border")
K_obs <- Kinhom(gorillas, lambda = pred.im, correction = "border")

# simulate the envelopes/bounds
K_env_ipp <- envelope(gorillas, fun = Kinhom, simulate = expression(rpoispp(lambda = pred_ipp.im)))
K_env <- envelope(gorillas, fun = Kinhom, simulate = expression(rpoispp(lambda = pred.im)))

# fit an IPP for comparison
m_ipp_poly <- gam(pt/wt ~ poly(elevation, 2) + poly(waterdist, 2) + poly(slopeangle, 2) + heat + slopetype + vegetation, data=dat, family=poisson(), weights=wt, method="REML")
# m_poly <- gam(pt/wt ~ poly(elevation, 2) + poly(waterdist, 2) + poly(slopeangle, 2) + heat + slopetype + vegetation + te(x, y, bs = "gp", k = 400, m = list(c(3, opt$minimum)), d = 2), data=dat, family=poisson(), weights=wt, method="REML")
m_poly <- gam(pt/wt ~ poly(elevation, 2) + poly(waterdist, 2) + poly(slopeangle, 2) + heat + slopetype + vegetation + s(x, y, bs = "gp", k = 400, m = c(3, opt$minimum)), data=dat, family=poisson(), weights=wt, method="REML")

# fit an IPP for comparison
m_ipp_sm <- gam(pt/wt ~ s(elevation) + s(waterdist) + s(slopeangle) + heat + slopetype + vegetation, data=dat, family=poisson(), weights=wt, method="REML")
# m_sm <- gam(pt/wt ~ s(elevation) + s(waterdist) + s(slopeangle) + heat + slopetype + vegetation + te(x, y, bs = "gp", k = 400, m = list(c(3, opt$minimum)), d = 2), data=dat, family=poisson(), weights=wt, method="REML")
m_sm <- gam(pt/wt ~ s(elevation) + s(waterdist) + s(slopeangle) + heat + slopetype + vegetation + s(x, y, bs = "gp", k = 400, m = c(3, opt$minimum)), data=dat, family=poisson(), weights=wt, method="REML")

# check the various information criteria
info_crit <- function(x){c(`GAM score`=ifelse(is.null(x$gcv.ubre),NA,x$gcv.ubre), logLik = logLik(x), AIC = AIC(x), BIC = AIC(x, k = log(nrow(pp))))}
data.frame(
  IPP_LINE = info_crit(m_ipp),
  IPP_POLY = info_crit(m_ipp_poly),
  IPP_SMOO = info_crit(m_ipp_sm),
  LGCP_LINE = info_crit(m),
  LGCP_POLY = info_crit(m_poly),
  LGCP_SMOO = info_crit(m_sm)
)

# Create Figure 2 in the manuscript:
plot.res <- 500
png(filename = paste0(home.wd, "/app_results_kfuncs.png"), res = plot.res, width = 6.2 * plot.res, height = 4.5 * plot.res)
layout(mat = matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE), widths = c(0.55, 0.45), heights = c(0.5, 0.5))
par(mar = c(2.1, 3.1, 2.1, 0))
plot(K_env_ipp$r, K_env_ipp$mmean, type = "n", ylim = range(c(K_env_ipp$obs, K_env_ipp$hi, K_env_ipp$lo)),
     ylab = "", xlab = "", xaxt = "n", yaxt = "n")
axis(side = 2, at = seq(0, 3e6, by = 1e6), labels = c("0", seq(1e6, 3e6, by = 1e6)))
polygon(c(rev(K_env_ipp$r), K_env_ipp$r), c(rev(K_env_ipp$hi), K_env_ipp$lo), col = 'grey80', border = NA)
lines(K_env_ipp$r, K_env_ipp$mmean, lty = "dashed")
lines(K_obs_ipp$r, K_obs_ipp$border, col = "red")
mtext(text = "A: Poisson Process", side = 3, cex = 1, line = 0.5, adj = 0)
par(mar = c(2.1, 0, 2.1, 0))
plot(log(pred_ipp.im), box = F, main = "", col = terrain.colors, ribbon = F)
points(pp[,c("x","y")], col = rgb(0,0,0,alpha = 0.25), pch = 4)
par(mar = c(3.1, 3.1, 1.1, 0))
plot(K_env$r, K_env$mmean, type = "n", ylim = range(c(K_env$obs, K_env$hi, K_env$lo)),
     ylab = "", xlab = "", yaxt = "n")
axis(side = 2, at = seq(0, 3e6, by = 1e6), labels = c("0", seq(1e6, 3e6, by = 1e6)))
polygon(c(rev(K_env$r), K_env$r), c(rev(K_env$hi), K_env$lo), col = 'grey80', border = NA)
lines(K_env$r, K_env$mmean, lty = "dashed")
lines(K_obs$r, K_obs$border, col = "red")
mtext(text = "distance (m)", side = 1, srt = 90, cex = 1, line = 2)
mtext(text = "Inhomogeneous K Function", side = 2, srt = 90, cex = 1, xpd = T, outer = T, line = -1)
mtext(text = "Estimated Intensity Surfaces with nest locations (x)", side = 4, srt = 90, cex = 1, xpd = T, outer = T, line = -16)
mtext(text = "B: log-Gaussian Cox Process", side = 3, cex = 1, line = 0.5, adj = 0)
legend(x = 0, y = 4e6, legend = c("Observed", "Theoretic", "95% Sim. Bounds"),
       col = c("red", "black", "grey80"), lty = c("solid", "dashed", "solid"), cex = 1,
       lwd = c(2, 2, 2), bty = "n")
par(mar = c(3.1, 0, 1.1, 0))
plot(log(pred.im), box = F, main = "", col = terrain.colors, ribbon = F)
points(pp[,c("x","y")], col = rgb(0,0,0,alpha = 0.25), pch = 4)
dev.off()

# demonstration using our gam LGCP wrapper function
source("gam_lgcp.R")

set.seed(1)
quad.id <- sample(1:nrow(quad), 10000)
set.seed(NULL)
subdat <- merge(pp, quad[quad.id, ], all = T)

time1 <- system.time(assign("mod", gam_lgcp(pt ~ elevation, data=subdat, weights=wt, range.interval = c(50,700))))
summary(mod)
time1
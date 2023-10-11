# make job array script

fit_model <- c("mgcv", "inla")
env_range <- c(20, 40, 60)
lat_range <- c(20, 40, 60)
sim <- 1:100

# expand out all combinations
tab <- data.frame(expand.grid(env_range, lat_range, sim, fit_model))
colnames(tab) <- c("env_range", "lat_range", "sim", "fit_model")
tab$job <- 1:nrow(tab)

write.csv(tab, file = "job_array.csv", row.names = F)

fit_model <- c("mgcv", "gam_opt", "youngman", "inla")
intercepts <- c(-4.5, -3.5, -2.5)
env_range <- c(10, 30, 50)
lat_range <- c(10, 30, 50)
sim <- 1:100

# expand out all combinations
tab <- data.frame(expand.grid(intercepts, env_range, lat_range, sim, fit_model))
colnames(tab) <- c("intercept", "env_range", "lat_range", "sim", "fit_model")
tab$job <- 1:nrow(tab)

write.csv(tab, file = "job_array_all.csv", row.names = F)

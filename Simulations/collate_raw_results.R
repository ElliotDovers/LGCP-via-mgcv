home.wd <- getwd()

# initialise the result storage
res <- NULL

# change into appropriate result folder
setwd(paste(home.wd, "Results", sep = "/"))
# get a list of all the individual result files
res.list <- list.files()[grepl("res_", list.files(), fixed = T)]
# inner loop through individual sim-by-model files
for (job in res.list) {
  # load in the individual simulation results
  load(job)
  # add to the data
  res <- rbind(res, res_tab)
  rm(res_tab)
}

setwd(home.wd)
save(list = "res", file = "collated_results.RDATA")
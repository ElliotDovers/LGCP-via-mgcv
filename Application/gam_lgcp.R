gam_lgcp <- function(formula, data, weights = NULL, coord.names = c("x", "y"), k = 100, range.interval, opt.tolerance = 3,
                     subset = NULL, na.action, offset = NULL, 
                     optimizer = c("outer", "newton"), control = list(), scale = 0, 
                     select = FALSE, knots = NULL, sp = NULL, min.sp = NULL, H = NULL, 
                     gamma = 1, fit = TRUE, paraPen = NULL, G = NULL, in.out = NULL, 
                     drop.unused.levels = TRUE, drop.intercept = NULL, nei = NULL, 
                     discrete = FALSE, ...) {
  
  mc <- match.call() # gets the arguments (must be updated for LGCP as below)
  call.list <- as.list(mc)
  
  # check the form of the weights
  object.supplied <- tryCatch(!is.null(weights), error = function(e) FALSE)
  if (!object.supplied) {
    weight.name <- deparse(substitute(weights))
    if (weight.name == "NULL") {
      stop("weights must be supplied for fitting a LGCP.\n\nThese must be quadrature weights in rows where the formula response = 0,\nThese must be some small number (e.g. 1e-6) in rows where the formula response = 1.")
    } else {
      wt.vec <- as.vector(data[,weight.name])
    }
  } else {
    wt.vec <- weights
  }
  # checks
  if ((!all(coord.names %in% colnames(data)))) {
    stop(paste0("One of 'coord.names', ", paste(coord.names, collapse = " or "), ", not found 'data'"))
  }

  # get the response variable out of the formula
  resp <- all.vars(formula[[2]])
  data$new.response <- data[ , resp] / wt.vec
  
  # alter the call according to requirements for a LGCP
  call.list$family <- poisson()
  call.list$data <- data
  call.list$weights <- wt.vec
  call.list$method <- "REML"
  # update the formula for an initial fit
  call.list$formula <- as.formula(paste0("new.response ~ ", as.character(formula)[3], " + s(", paste(coord.names, collapse = ", "), ", bs=\"gp\", k=", deparse(k), ", m=3)"))
  # remove the function of the call
  call.list[[1]] <- NULL
  # fit an initial model to obtain warm starting parameters
  init.mod <- do.call(mgcv::gam, call.list)
  warm.starts <- init.mod$coefficients
  warm.starts[(length(init.mod$coefficients) - k + 1):length(init.mod$coefficients)] <- 0
  # update the starting parameters
  call.list$start <- warm.starts
  # set up the object function to be optimized
  objective_fn <- function(rho) {
    call.list$formula <- as.formula(paste0("new.response ~ ", as.character(formula)[3], " + s(", paste(coord.names, collapse = ", "), ", bs=\"gp\", k=", deparse(k), ", m=c(3,", deparse(rho), "))"))
    tmp.m = do.call(mgcv::gam, call.list)
    return(tmp.m$gcv.ubre) # the "method" specific criterion
  }
  # calculate the optimum
  opt <- optimize(objective_fn, interval = range.interval, tol = opt.tolerance)
  # adjust the formula for the optimized spatial range parameter
  call.list$formula <- as.formula(paste0("new.response ~ ", as.character(formula)[3], " + s(", paste(coord.names, collapse = ", "), ", bs=\"gp\", k=", deparse(k), ", m=c(3,", deparse(opt$minimum), "))"))
  # fit the final model
  res <- do.call(mgcv::gam, call.list)
  return(res)
}

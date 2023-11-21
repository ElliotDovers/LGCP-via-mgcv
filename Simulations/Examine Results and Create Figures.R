## Analyse Scenario ##

home.wd <- getwd()

# load up the pre-collated results object
load("collated_results.RDATA")

library(dplyr)
library(ggplot2)
library(gridExtra)
library(xtable)

plot.res <- 500

# define the scenarios
res$scen <- paste(res$env_range, res$lat_range, sep = ":")
res$scenario <- factor(res$scen,
                       levels = c("10:10", "10:30", "10:50", "30:10", "30:30", "30:50", "50:10", "50:30", "50:50"),
                       labels = c("rho[X]==10~~~~rho[xi]==10",
                                  "rho[X]==10~~~~rho[xi]==30",
                                  "rho[X]==10~~~~rho[xi]==50",
                                  "rho[X]==30~~~~rho[xi]==10",
                                  "rho[X]==30~~~~rho[xi]==30",
                                  "rho[X]==30~~~~rho[xi]==50",
                                  "rho[X]==50~~~~rho[xi]==10",
                                  "rho[X]==50~~~~rho[xi]==30",
                                  "rho[X]==50~~~~rho[xi]==50")
)
res$scenario_lat <- factor(res$lat_range,
                           levels = c("10", "30", "50"),
                           labels = c("rho==10",
                                      "rho==30",
                                      "rho==50")
)
res$scenario_env <- factor(res$env_range,
                           levels = c("10", "30", "50"),
                           labels = c("rho[X]==10",
                                      "rho[X]==30",
                                      "rho[X]==50")
)

# OLD SET OF APPROACHES TRIALLED MORE COMBINATIONS OF CRITERIA - UPDATED FOR S() SIMULATIONS
# res$Approach <- factor(paste(res$GP_APPROX, res$CRIT, res$method, sep = ":"),
#                     levels = c('gp.smooth:NA:REML', 'gp.smooth:ll:REML', 'gp.smooth:aic:REML', 'gp.smooth:score:REML', 'gp.smooth:bic:REML', 'gp.smooth:NA:GCV.Cp', 'gp.smooth:ll:GCV.Cp', 'gp.smooth:aic:GCV.Cp', 'gp.smooth:score:GCV.Cp', 'gp.smooth:bic:GCV.Cp', 'tprs.smooth:NA:REML', 'tprs.smooth:NA:GCV.Cp', 'pcmatern:NA:NA', 'default:NA:NA'),
#                     labels = c("mgcv REML GP Def.", "mgcv REML GP LL", "mgcv REML GP AIC", "mgcv REML GP SCORE", "mgcv REML GP BIC",  "mgcv GCV GP Def.", "mgcv GCV GP LL", "mgcv GCV GP AIC", "mgcv GCV GP SCORE", "mgcv GCV GP BIC", "mgcv REML TPRS Def.", "mgcv GCV TPRS Def.",  "INLA Fuglstad P.C.", "INLA Def.")
# )
res$Approach <- factor(paste(res$GP_APPROX, res$CRIT, res$method, sep = ":"),
                       levels = c('gp.smooth:NA:REML', 'gp.smooth:score:REML', 'gp.smooth:NA:GCV.Cp', 'gp.smooth:score:GCV.Cp', 'tprs.smooth:NA:REML', 'tprs.smooth:NA:GCV.Cp', 'pcmatern:NA:NA'),
                       labels = c("mgcv REML GP Def.", "mgcv REML GP SCORE", "mgcv GCV GP Def.", "mgcv GCV GP SCORE", "mgcv REML TPRS Def.", "mgcv GCV TPRS Def.",  "INLA Fuglstad P.C.")
)

# put together results for appendices
tab <- res %>% filter(scen == "30:30" & intercept == -3.5) %>% group_by(Approach, K) %>% summarise(MAE = mean(MAE), RMSE_BETA = sqrt(mean(SQER_BETA)), COVERAGE = mean(COVER_BETA), TIMING = mean(ALL_TIME), RMSE_RHO = sqrt(mean(SQER_RHO)), NSIMS = length(sim))
# remove some of the unnecessary comparisons
tab <- tab %>% filter(!Approach %in% c("mgcv REML GP BIC", "mgcv GCV GP BIC", "mgcv REML TPRS Def."))
# tidy up some columns
tab$K <- factor(tab$K)
tab$Software <- "mgcv"
tab$Software[tab$Approach %in% c("INLA Fuglstad P.C.", "INLA Def.")] <- "INLA"
tab$Method <- "REML"
tab$Method[grepl("GCV", tab$Approach)] <- "UBRE"
tab$Method[tab$Approach == "INLA Fuglstad P.C."] <- "P.C. Priors"
tab$Method[tab$Approach == "INLA Def."] <- "Default"
tab$Basis_Functions <- ""
tab$Basis_Functions[grepl(" GP ", tab$Approach)] <- "GP"
tab$Basis_Functions[grepl(" TPRS ", tab$Approach)] <- "TPRS"
tab$Basis_Functions[tab$Approach %in% c("INLA Fuglstad P.C.", "INLA Def.")] <- "SPDE Mesh"
tab$Range_Method <- ""
tab$Range_Method[grepl(" AIC", tab$Approach)] <- "AIC"
tab$Range_Method[grepl(" LL", tab$Approach)] <- "log-Lik."
tab$Range_Method[grepl(" SCORE", tab$Approach)] <- "REML/UBRE"
tab$Range_Method[tab$Approach %in% c("INLA Fuglstad P.C.", "INLA Def.")] <- "Default"
tab$COVERAGE <- paste0(tab$COVERAGE * 100, "%")
tab$RMSE_RHO[tab$Range_Method == ""] <- NA
print(
  xtable(
    tab[, c("Software", "Method", "Basis_Functions", "Range_Method", "K", "MAE", "COVERAGE", "TIMING", "RMSE_BETA", "RMSE_RHO")], label = "append:tab:comp_sim", digits = 3, caption.placement = "top", caption = "Result comparisons for models used to fit LGCPs."
  ) , include.rownames = F
)

# turn K into a factor
res$k <- factor(res$K,
                levels = c("25", "100", "200", "300", "400"),
                labels = c("25 (mgcv Default)", "100", "200", "300", "400")
)
res$k <- res$K
res$k[!is.na(res$INLA_K)] <- res$INLA_K[!is.na(res$INLA_K)]
# create a factor for intercept that is expected N
res$E_N <- factor(res$intercept,
                  levels = sort(unique(res$intercept)),
                  labels = paste0("E(N)==", round((res %>% group_by(intercept) %>% summarise(E_N = mean(n)))$E_N))
)

# set up colors for all models - UPDATED FOR S() SIMULATIONS
# fit_cols_all <- c("darkkhaki", "darkolivegreen1", "darkolivegreen2", "darkolivegreen3", "darkolivegreen4", # Our mgcv models REML
#                   "darkgoldenrod1", "khaki", "darkgoldenrod2", "darkgoldenrod4", "purple", # Our mgcv models GCV (ACUTALLY UBRE)
#                   "tomato1", "tomato4", # Youngman et all model
#                   "royalblue1", "royalblue4" # INLA models
# )
fit_cols_all <- c("darkkhaki", "darkolivegreen3", "darkgoldenrod1", "darkgoldenrod4", "tomato1",
                  "tomato4", "royalblue1"
)
# separate out the results for the INLA comparison
# 1: "mgcv REML GP Def."
# 2: "mgcv REML GP LL"
# 3: "mgcv REML GP AIC"
# 4: "mgcv REML GP SCORE"
# 5: "mgcv REML GP BIC"
# 6: "mgcv GCV GP Def."
# 7: "mgcv GCV GP LL"
# 8: "mgcv GCV GP AIC"
# 9: "mgcv GCV GP SCORE"
# 10: "mgcv GCV GP BIC"
# 11: "mgcv REML TPRS Def."
# 12: "mgcv GCV TPRS Def."
# 13: "INLA Fuglstad P.C."
# 14: "INLA Def."

### THESE HAVE CHANGED WITH THE REDUCED CRITERIA COMPARISONS (UPDATED FOR S() SIMULATIONS):
# 1: "mgcv REML GP Def." (was 1)
# 2: "mgcv REML GP SCORE" (was 4)
# 3: "mgcv GCV GP Def." (was 6)
# 4: "mgcv GCV GP SCORE" (was 9)
# 5: "mgcv REML TPRS Def." (was 11)
# 6: "mgcv GCV TPRS Def." (was 12)
# 7: "INLA Fuglstad P.C." (was 13)

# mods_to_comp <- c(6,13) UPDATED FOR S() SIMULATIONS
mods_to_comp <- c(3,7)
res_sec1 <- res %>% filter(env_range == 30 & lat_range == 30 & intercept == -3.5 & Approach %in% levels(res$Approach)[mods_to_comp])

# define the model colours for plotting
fit_cols <- fit_cols_all[mods_to_comp]

# transform the results to show means and SD
res_sec1 %>% group_by(Approach, K) %>% summarise(y = mean(MAE), sd = sd(MAE))

p1 <- ggplot(data = res_sec1 %>% group_by(Approach, K, k) %>% summarise(y = mean(MAE), sd = sd(MAE)), aes(x = K, y = y, fill = Approach, color = Approach)) +
  geom_point(position = position_dodge(3)) + geom_errorbar(aes(ymin = y-sd, ymax = y + sd), position = position_dodge(3), width = 15) + geom_line(position = position_dodge(3)) +
  scale_y_continuous(name = expression(paste(MAE,": |", lambda, " - ", hat(lambda), "|"))) +
  scale_x_continuous(name = "") +
  scale_color_manual(values = fit_cols) +
  scale_fill_manual(values = alpha(fit_cols, alpha = 0.25)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), legend.position = "none",
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 67, b = 0, l = 0)
  ) +
  ggtitle(label = "A")

cps <- res_sec1 %>% group_by(Approach, K) %>% summarise(cp = mean(COVER_BETA)*100)
p2 <- ggplot(data = cps, aes(x = K, y = cp, color = Approach, fill = Approach)) +
  geom_point(position = position_dodge(3))+ geom_line(position = position_dodge(3)) +
  geom_abline(slope = 0, intercept = log(95), col = "red", lty = "dashed") +
  # geom_abline(slope = 0, intercept = 95, col = "red", lty = "dashed") +
  scale_y_continuous(name = expression(paste("Cover. Prob. ", beta[1], " (", alpha, " = 0.05)")), trans = "log", breaks = c(50, 60, 70, 80, 90, 100), labels = c("50%", "60%", "70%", "80%", "90%", "100%")) +
  # scale_y_continuous(name = expression(paste("Cover. Prob. ", hat(beta)[1], " (", alpha, " = 0.05)")), breaks = c(50, 60, 70, 80, 90, 95, 100), labels = c("50%", "60%", "70%", "80%", "90%", "95%", "100%")) +
  scale_x_continuous(name = "") +
  scale_color_manual(values = fit_cols, name = "Software", labels = c("mgcv", "INLA")) +
  scale_fill_manual(values = alpha(fit_cols, alpha = 0.5), name = "Software", labels = c("mgcv", "INLA")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.title.y = element_text(vjust = 0),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.key = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 3)
  ) +
  ggtitle(label = "B")

p3 <- ggplot(data = res_sec1 %>% group_by(Approach, K) %>% summarise(y = mean(ALL_TIME), sd = sd(ALL_TIME)), aes(x = K, y = y, color = Approach, fill = Approach)) +
  geom_point(position = position_dodge(3)) + geom_errorbar(aes(ymin = y-sd, ymax = y + sd), position = position_dodge(3), width = 15) + geom_line(position = position_dodge(3)) +
  scale_y_continuous(name = "Comp. Time (sec)") +
  scale_x_continuous(name = "Basis Dimension (mgcv knots/INLA vertices)", breaks = c(25, 100, 200, 300, 400), labels = c("25 (mgcv def.)", "100", "200", "300", "400")) +
  scale_color_manual(values = fit_cols, name = "Software", labels = c("INLA", "mgcv")) +
  scale_fill_manual(values = alpha(fit_cols, alpha = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.title.y = element_text(vjust = 7),
        plot.margin = margin(t = 0, r = 67, b = 6, l = 20), legend.position = "none"
  ) +
  ggtitle(label = "C")

# Main result plot
png(filename = paste0(home.wd, "/Figures/mgcv_inla_sim_results.png"), width = 6.2 * plot.res, height = 5.5 * plot.res, res = plot.res)
grid.arrange(p1,p2,p3, nrow = 3)
dev.off()

## Appendices Plots ##

# full results plots
# mods_to_comp <- c(6, 4, 9, 13) UPDATED FOR S() SIMULATIONS
mods_to_comp <- c(3, 2, 4, 7)
tmp.plot.data <- res %>% filter(env_range == 30 & Approach %in% levels(res$Approach)[mods_to_comp]) %>% group_by(scenario_lat, E_N, Approach, K) %>% summarise(y = mean(MAE), sd = sd(MAE))
tmp.plot.data$Approach <- factor(tmp.plot.data$Approach, levels = levels(res$Approach)[mods_to_comp])
png(filename = paste0(home.wd, "/Figures/append_sim_results_mae.png"), width = 6.3 * plot.res, height = 5.5 * plot.res, res = plot.res)
ggplot(data = tmp.plot.data, aes(x = K, y = y, color = Approach, fill = Approach)) +
  # geom_point(position = position_dodge(5)) + geom_errorbar(aes(ymin = y-sd, ymax = y + sd), position = position_dodge(5), width = 30) + geom_line(position = position_dodge(5)) +
  geom_point(position = position_dodge(3)) + geom_line(position = position_dodge(3)) +
  facet_wrap( ~ scenario_lat + E_N, labeller = label_parsed, scales = "free_y") +
  scale_y_continuous(name = expression(paste(MAE,": |", lambda, " - ", hat(lambda), "|"))) +
  scale_x_continuous(name = "Basis Dimension (mgcv knots/INLA vertices)", breaks = c(25, 100, 200, 300, 400), labels = c("25", "100", "200", "300", "400")) +
  scale_color_manual(values = fit_cols_all[mods_to_comp], name = "Approach", labels = c("mgcv Default", expression(paste("mgcv ", hat(rho), " via REML")), expression(paste("mgcv ", hat(rho), " via UBRE")), "INLA")) +
  scale_fill_manual(values = alpha(fit_cols_all[mods_to_comp], alpha = 0.5), name = "Approach", labels = c("mgcv Default", expression(paste("mgcv ", hat(rho), " via REML")), expression(paste("mgcv ", hat(rho), " via UBRE")), "INLA")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.title.y = element_text(vjust = 0), legend.text.align = 0,
        # axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )
dev.off()

tmp.plot.data <- res %>% filter(env_range == 30 & Approach %in% levels(res$Approach)[mods_to_comp]) %>% group_by(scenario_lat, E_N, Approach, K) %>% summarise(y = mean(COVER_BETA))
tmp.plot.data$Approach <- factor(tmp.plot.data$Approach, levels = levels(res$Approach)[mods_to_comp])
png(filename = paste0(home.wd, "/Figures/append_sim_results_cover.png"), width = 6.2 * plot.res, height = 5.5 * plot.res, res = plot.res)
ggplot(data = tmp.plot.data, aes(x = K, y = y, color = Approach, fill = Approach)) +
  geom_point(position = position_dodge(3)) + geom_line(position = position_dodge(3)) +
  facet_wrap( ~ scenario_lat + E_N, labeller = label_parsed)+#, scales = "free_y") +
  scale_y_continuous(name = expression(paste("Cover. Prob. ", beta[1], " (", alpha, " = 0.05)"))) +
  scale_x_continuous(name = "Basis Dimension (mgcv knots/INLA vertices)", breaks = c(25, 100, 200, 300, 400), labels = c("25", "100", "200", "300", "400")) +
  scale_color_manual(values = fit_cols_all[mods_to_comp], name = "Approach", labels = c("mgcv Default", expression(paste("mgcv ", hat(rho), " via REML")), expression(paste("mgcv ", hat(rho), " via UBRE")), "INLA")) +
  scale_fill_manual(values = alpha(fit_cols_all[mods_to_comp], alpha = 0.5), name = "Approach", labels = c("mgcv Default", expression(paste("mgcv ", hat(rho), " via REML")), expression(paste("mgcv ", hat(rho), " via UBRE")), "INLA")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.title.y = element_text(vjust = 0), legend.text.align = 0,
        # axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
  ) +
  geom_abline(slope = 0, intercept = 0.95, col = "black", lty = "dashed")
dev.off()

tmp.plot.data <- res %>% filter(env_range == 30 & Approach %in% levels(res$Approach)[mods_to_comp]) %>% group_by(scenario_lat, E_N, Approach, K) %>% summarise(y = mean(TIME), sd = sd(TIME))
tmp.plot.data$Approach <- factor(tmp.plot.data$Approach, levels = levels(res$Approach)[mods_to_comp])
png(filename = paste0(home.wd, "/Figures/append_sim_results_comptime.png"), width = 6.2 * plot.res, height = 5.5 * plot.res, res = plot.res)
ggplot(data = tmp.plot.data, aes(x = K, y = y, color = Approach, fill = Approach)) +
  geom_point(position = position_dodge(3)) + geom_line(position = position_dodge(3)) +
  facet_wrap( ~ scenario_lat + E_N, labeller = label_parsed)+#, scales = "free_y") +
  scale_y_continuous(name = "Comp. Time (sec)", trans = "log", breaks = c(5, 15, 30, 60, 120, 240)) +
  scale_x_continuous(name = "Basis Dimension (mgcv knots/INLA vertices)", breaks = c(25, 100, 200, 300, 400), labels = c("25", "100", "200", "300", "400")) +
  scale_color_manual(values = fit_cols_all[mods_to_comp], name = "Approach", labels = c("mgcv Default", expression(paste("mgcv ", hat(rho), " via REML")), expression(paste("mgcv ", hat(rho), " via UBRE")), "INLA")) +
  scale_fill_manual(values = alpha(fit_cols_all[mods_to_comp], alpha = 0.5), name = "Approach", labels = c("mgcv Default", expression(paste("mgcv ", hat(rho), " via REML")), expression(paste("mgcv ", hat(rho), " via UBRE")), "INLA")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.title.y = element_text(vjust = 0), legend.text.align = 0,
        # axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )
dev.off()
        

# Choosing k plot for appendices
# mods_to_comp <- c(6,4,9) UPDATED FOR S() SIMULATIONS
mods_to_comp <- c(3, 2, 4)
tmp.plot.data <- res %>% filter(env_range == 30 & Approach %in% levels(res$Approach)[mods_to_comp]) %>% group_by(scenario_lat, E_N, Approach, K) %>% summarise(y = mean(EDF), sd = sd(EDF))
tmp.plot.data$Approach <- factor(tmp.plot.data$Approach, levels = levels(res$Approach)[mods_to_comp])
png(filename = paste0(home.wd, "/Figures/append_sim_results_choosing_k.png"), width = 6.2 * plot.res, height = 5.5 * plot.res, res = plot.res)
ggplot(data = tmp.plot.data, aes(x = K, y = y, color = Approach, fill = Approach, linetype = "dummy")) +
  geom_point(position = position_dodge(15)) + geom_errorbar(aes(ymin = y-sd, ymax = y + sd, linetype = NULL), position = position_dodge(15), width = 30) +
  facet_wrap( ~ scenario_lat + E_N, labeller = label_parsed, scales = "free_y") +
  scale_y_continuous(name = "Effective Degrees of Freedom") +
  scale_x_continuous(name = "Basis Dimension (mgcv knots/INLA vertices)", breaks = c(25, 100, 200, 300, 400), labels = c("25", "100", "200", "300", "400")) +
  scale_color_manual(values = fit_cols_all[mods_to_comp], name = "Estimated Range", labels = c("Default", "Opt. via REML", "Opt. via UBRE")) +
  scale_fill_manual(values = alpha(fit_cols_all[mods_to_comp], alpha = 0.5), name = "Estimated Range", labels = c("Default", "Opt. via REML", "Opt. via UBRE")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.title.y = element_text(vjust = 0),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0), #legend.position = "none"
  ) +
  geom_abline(slope = 0.25, intercept = 0, col = "black", linetype = "dashed") +
  scale_linetype_manual(name = expression(paste("Ratio: ", frac(k, edf))), values = c("dashed"), labels = "4x")
dev.off()

# Estimating spatial range parameter plot for appendices
# mods_to_comp <- c(4, 9, 13) UPDATED FOR S() SIMULATIONS
mods_to_comp <- c(2, 4, 7)
true_ranges <- expand.grid(levels(res$scenario_lat), levels(res$E_N))
colnames(true_ranges) <- c("scenario_lat", "E_N")
true_ranges$intercepts <- rep(c(10, 30, 50), times = 3)
true_ranges$type <- "rho"
tmp.plot.data <- res %>% filter(env_range == 30 & Approach %in% levels(res$Approach)[mods_to_comp]) %>% group_by(scenario_lat, E_N, Approach, K) %>% summarise(y = mean(RHO_HAT), sd = sd(RHO_HAT))
tmp.plot.data$Approach <- factor(tmp.plot.data$Approach, levels = levels(res$Approach)[mods_to_comp])
png(filename = paste0(home.wd, "/Figures/append_sim_results_estimated_rho.png"), width = 6.2 * plot.res, height = 5.5 * plot.res, res = plot.res)
ggplot(data = tmp.plot.data, aes(x = K, y = y, color = Approach, fill = Approach)) +
  geom_point(position = position_dodge(15)) + geom_errorbar(aes(ymin = y-sd, ymax = y + sd), position = position_dodge(15), width = 30) +
  facet_wrap( ~ scenario_lat + E_N, labeller = label_parsed, scales = "free_y") +
  scale_y_continuous(name = expression(hat(rho)[xi])) +
  scale_x_continuous(name = "Basis Dimension (mgcv knots/INLA vertices)", breaks = c(25, 100, 200, 300, 400), labels = c("25", "100", "200", "300", "400")) +
  scale_color_manual(values = fit_cols_all[mods_to_comp], name = "Approach", labels = c("mgcv via REML", "mgcv via UBRE", "INLA")) +
  scale_fill_manual(values = alpha(fit_cols_all[mods_to_comp], alpha = 0.5), name = "Approach", labels = c("mgcv via REML", "mgcv via UBRE", "INLA")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.title.y = element_text(vjust = 0),
        # axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
  ) +
  geom_hline(
    data=true_ranges,
    mapping=aes(yintercept=intercepts, linetype = type),
    size=0.5, color="black",# lty = "dashed",
    key_glyph="path") +
  scale_linetype_manual(name = "Simulated\nSpatial\nRange", values = c("dashed"), labels = expression(paste("  ", rho)))
dev.off()

# # REML vs UBRE
# mods_to_comp <- c(1, 3)
# tmp.plot.data <- res %>% filter(env_range == 30 & Approach %in% levels(res$Approach)[mods_to_comp]) %>% group_by(scenario_lat, E_N, Approach, K) %>% summarise(y = mean(COVER_BETA), sd = sd(COVER_BETA))
# tmp.plot.data$Approach <- factor(tmp.plot.data$Approach, levels = levels(res$Approach)[mods_to_comp])
# png(filename = paste0(home.wd, "/Figures/append_sim_results_reml_vs_ubre.png"), width = 6.2 * plot.res, height = 5.5 * plot.res, res = plot.res)
# ggplot(data = tmp.plot.data, aes(x = K, y = y, color = Approach, fill = Approach)) +
#   geom_point(position = position_dodge(75)) + geom_line(position = position_dodge(75)) +
#   facet_wrap( ~ scenario_lat + E_N, labeller = label_parsed)+#, scales = "free_y") +
#   scale_y_continuous(name = expression(paste("Cover. Prob. ", hat(beta)[1], " (", alpha, " = 0.05)"))) +
#   scale_x_continuous(name = "# Basis Functions (mgcv knots/INLA vertices)", breaks = c(25, 100, 200, 300, 400), labels = c("25", "100", "200", "300", "400")) +
#   scale_color_manual(values = fit_cols_all[mods_to_comp], name = "Method", labels = c("mgcv Default", expression(paste("mgcv Opt. ", hat(rho)[xi], " via REML")), expression(paste("mgcv Opt. ", hat(rho)[xi], " via UBRE")), "INLA")) +
#   scale_fill_manual(values = alpha(fit_cols_all[mods_to_comp], alpha = 0.5), name = "Approach", labels = c("mgcv Default", expression(paste("mgcv Opt. ", hat(rho)[xi], " via REML")), expression(paste("mgcv Opt. ", hat(rho)[xi], " via UBRE")), "INLA")) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
#         axis.text.y = element_text(size=10), axis.title.y = element_text(vjust = 0), legend.text.align = 0,
#         # axis.text.x = element_blank(), axis.ticks.x = element_blank(),
#         plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
#   ) +
#   geom_abline(slope = 0, intercept = 0.95, col = "black", lty = "dashed")
# dev.off()

# GP vs TPRS
# mods_to_comp <- c(1,6,11,12) UPDATED FOR S() SIMULATIONS
mods_to_comp <- c(1, 3, 5, 6)
res_append <- res %>% filter(env_range == 30 & lat_range == 30 & intercept == -3.5 & method == "REML" & !is.na(method) & Approach %in% levels(res$Approach)[mods_to_comp])
fit_cols <- c("darkgoldenrod1", "purple")

p1 <- ggplot(data = res_append %>% group_by(Approach, K) %>% summarise(y = mean(MAE), sd = sd(MAE)), aes(x = K, y = y, fill = Approach, color = Approach)) +
  geom_point(position = position_dodge(3)) + geom_errorbar(aes(ymin = y-sd, ymax = y + sd), position = position_dodge(3), width = 15) + geom_line(position = position_dodge(3)) +
  scale_y_continuous(name = expression(paste(MAE,": |", lambda, " - ", hat(lambda), "|"))) +
  scale_x_continuous(name = "") +
  scale_color_manual(values = fit_cols) +
  scale_fill_manual(values = alpha(fit_cols, alpha = 0.25)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), legend.position = "none",
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 100, b = 0, l = 0)
  ) +
  ggtitle(label = "A")


cps <- res_append %>% group_by(Approach, K) %>% summarise(cp = mean(COVER_BETA)*100)
p2 <- ggplot(data = cps, aes(x = K, y = cp, color = Approach)) +
  geom_point()+ geom_line() +
  geom_abline(slope = 0, intercept = log(95), col = "red", lty = "dashed") +
  # geom_abline(slope = 0, intercept = 95, col = "red", lty = "dashed") +
  scale_y_continuous(name = expression(paste("Cover. Prob. ", beta[1], " (", alpha, " = 0.05)")), trans = "log", breaks = c(50, 60, 70, 80, 89, 92, 95, 98, 100), labels = c("50%", "60%", "70%", "80%", "89%", "92%", "95%", "98%", "100%")) +
  # scale_y_continuous(name = expression(paste("Cover. Prob. ", hat(beta)[1], " (", alpha, " = 0.05)")), breaks = c(50, 60, 70, 80, 90, 95, 100), labels = c("50%", "60%", "70%", "80%", "90%", "95%", "100%")) +
  scale_x_continuous(name = "") +
  scale_color_manual(values = fit_cols, name = "Basis Functions", labels = c("GP", "TPRS")) +
  scale_fill_manual(values = alpha(fit_cols, alpha = 0.5), name = "Basis Functions", labels = c("GP", "TPRS")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.title.y = element_text(vjust = 1),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.key = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 3)
  ) +
  ggtitle(label = "B")

p3 <- ggplot(data = res_append %>% group_by(Approach, K) %>% summarise(y = mean(ALL_TIME), sd = sd(ALL_TIME)), aes(x = K, y = y, color = Approach, fill = Approach)) +
  geom_point(position = position_dodge(3)) + geom_errorbar(aes(ymin = y-sd, ymax = y + sd), position = position_dodge(3), width = 15) + geom_line(position = position_dodge(3)) +
  scale_y_continuous(name = "Comp. Time (sec)") +
  scale_x_continuous(name = "Basis Dimension (mgcv knots/INLA vertices)", breaks = c(25, 100, 200, 300, 400), labels = c("25 (mgcv def.)", "100", "200", "300", "400")) +
  scale_color_manual(values = fit_cols, name = "Software", labels = c("INLA", "mgcv")) +
  scale_fill_manual(values = alpha(fit_cols, alpha = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.title.y = element_text(vjust = 6),
        plot.margin = margin(t = 0, r = 100, b = 6, l = 15), legend.position = "none"
  ) +
  ggtitle(label = "C")

# Plot comparing results for TPRS vs GP basis functions
png(filename = paste0(home.wd, "/Figures/append_sim_results_gp_tprs.png"), width = 6.2 * plot.res, height = 5.5 * plot.res, res = plot.res)
grid.arrange(p1,p2,p3, nrow = 3)
dev.off()

################################################################################
# =============================================================================
# Motivation example introduction
# =============================================================================
devtools::load_all()

set.seed(2)
gen_dat = function(n) {
  x1 = runif(n)
  x2 = runif(n, min = 0, max = 1 - x1)
  data.frame(x1 = x1, x2 = x2)
}

df = gen_dat(100)

y.fun = function(X.model = NULL, newdata) {
  exp(newdata$x1 + newdata$x2)
}

y = y.fun(newdata = df)
pred = Predictor$new(predict.fun = y.fun, data = df, y = y)
fi = FeatureImp$new(pred, loss = "mse", compare = "difference")

grid.dat = expand.grid(x1 = seq(from = min(df$x1), to = max(df$x1), length.out = 20),
                       x2 = seq(from = min(df$x2), to = max(df$x2), length.out = 20))
grid.dat$predicted = y.fun(mod, grid.dat)

p = ggplot(df) +
  geom_tile(data = grid.dat, aes(x = x1, y = x2, fill = predicted)) +
  geom_point(aes(x = x1, y = x2), size = 3) +
  scale_fill_continuous("Prediction") +
  scale_x_continuous(TeX("Feature x_1")) +
  scale_y_continuous(TeX("Feature x_2")) +
  scale_fill_viridis(option = "D") +
  theme(legend.position = "top")
intercept = 0.9

max_pred = max(pred$predict(df)[[1]])
p1 = FeatureEffect$new(pred, "x1", method = "pdp")$plot(rug=FALSE) +
  #ggtitle(expression("PDP for"~x[1])) + xlab(expression(x[1])) +
  theme(plot.title = element_text(size = 16)) +
#  geom_hline(yintercept = max_pred, lty = 2) +
  geom_point(data = df, aes(x = x1, y = y), color = "darkgrey") +
  scale_x_continuous(TeX("Feature x_1")) +
	scale_y_continuous("Model prediction")

pdf(file = sprintf("%s/correlation-problem.pdf", fig_dir), width = 8, height = 4)
print(p + p1)
dev.off()


# =============================================================================
# M- Plot FAil
# =============================================================================

set.seed(42)
n = 500
x1 = rnorm(n, mean = 0, sd = 1)
x2 = x1 + rnorm(n, mean = 0, sd = 0.20)
x3 = rnorm(n, mean = 0, sd = 1)
x = data.frame(x1, x2, x3)

fx = function(mod=NULL, newdata){
  newdata$x1  - 0.1 * newdata$x2 + newdata$x3
}
predx = Predictor$new(predict.fun = fx, data = x, y = fx(newdata = x))
xgrid = seq(from = min(x$x1), to = max(x$x1), length.out = 30)
window_width = (max(x$x1) - min(x$x1)) / 20
mplot_data = mplot(predx, "x2", xgrid, window_width)$results
pdp_data = FeatureEffect$new(predx, "x2", method = "pdp")$results
plot_data = data.frame(rbind(pdp_data, mplot_data))
plot_data$.type[plot_data$.type == "pdp"] = "PDP"
plot_data$.type[plot_data$.type == "m-plot"] = "cond. PDP"

# Feature Importance
xPFI = FeatureImp$new(predx, loss = "mae", compare = "difference")$results
xPFI$type = "marginal"
results = data.frame(xPFI)

conds = fit_conditionals(predx$data$X)
xPFI =  grouped_pfi(predx, loss = Metrics::mae, conds = conds, repetitions = 5)
xPFI = data.frame(xPFI)
xPFI$type = "conditional"
cols = c("feature", "importance", "type")
results = rbind(results[,cols], xPFI[,cols])

# Combine plot and table
p1 = ggplot(plot_data) +
  geom_line(aes(x = x2, y = .value, group = .type, color = .type, lty = .type), size = 2) +
  scale_y_continuous("Average prediction") +
  scale_x_continuous(TeX("Feature x_2")) +
  scale_color_viridis("", discrete = TRUE, option = "D") +
  scale_linetype_discrete("") +
  theme(legend.position = "top")

results$type = factor(results$type)
results$feature = factor(results$feature, levels = c("x3", "x2", "x1"))
p2 = ggplot(results) +
  geom_col(aes(x = feature, y = importance, fill = type), position = position_dodge(1)) +
  coord_flip() +
  scale_x_discrete("") +
  scale_y_continuous("Importance") +
  scale_fill_viridis("", discrete = TRUE, option = "D") +
  theme(legend.position = "top")

pdf(file = sprintf("%s/mplot-fail.pdf", fig_dir), height = 4.2, width = 8)
print(p1 + p2)
dev.off()


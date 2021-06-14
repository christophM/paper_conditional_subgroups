# =============================================================================
# Motivation example introduction
# =============================================================================
devtools::load_all()

set.seed(1)
gen_dat = function(n) {
  x1 = rnorm(n)
  x2 = rnorm(n, mean = x1)
  data.frame(x1 = x1, x2 = x2)
}

x = gen_dat(5000)
x = round(x, 2)


# cor(x$x1, x$x2)
#0.7165793

y.fun = function(X.model = NULL, newdata) {
  newdata$x1
}

x$y = y.fun(newdata = x) + rnorm(nrow(x))

mod = lm(y ~ x1, data = x)

pred = Predictor$new(mod, data = x, y = "y")
xgrid = seq(from = min(x$x1), to = max(x$x1), length.out = 30)
window_width = (max(x$x1) - min(x$x1)) / 10
mplot_data = mplot(pred, "x2", xgrid, window_width)$results
pdp_data = FeatureEffect$new(pred, "x2", method = "pdp")$results
plot_data = data.frame(rbind(pdp_data, mplot_data))
plot_data$.type[plot_data$.type == "pdp"] = "PDP"
plot_data$.type[plot_data$.type == "m-plot"] = "cond. PDP"


yscale = scale_y_continuous("Prediction", limits = c(-2, 2))
# Combine plot and table
p1 = ggplot(plot_data) +
  geom_line(aes(x = x2, y = .value, group = .type, color = .type, lty = .type), size = 2) +
  yscale +
  scale_x_continuous(TeX("Feature x_2")) +
  scale_color_viridis("", discrete = TRUE, option = "D") +
  scale_linetype_discrete("") +
  theme(legend.position = "top")

# TODO: Add cs-PDPs

ctr = rpart::rpart.control(maxdepth = 2, minbucket = 30)
conds = fit_conditionals(pred$data$X, type = "cart", ctrl = ctr)
p2 = grouped_pdp(pred, "x2", conds[["x2"]])$p +
  scale_x_continuous(TeX("Feature x_2")) +
  scale_color_discrete("subgroups") +
  yscale
p = (p1 + p2)

ggsave(p, file = sprintf("%s/mplot-fail-again.pdf", fig_dir), height = 4, width = 8)

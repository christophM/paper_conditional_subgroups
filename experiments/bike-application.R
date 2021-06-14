# =============================================================================
# Bike rental application
# =============================================================================

library("mlr3")
library("mlr3learners")
library("here")
library("ranger")
library("parallel")
library("partykit")
library("rpart")

set.seed(2)
test_size = 0.3
devtools::load_all()

#==============================================================================
# Preparing the data
#==============================================================================
# Measure: Correlation, predictive dependence
# Load bike dataset
bike = read.csv(paste0(data_dir,"bike-sharing-daily.csv"))

bike$weekday = factor(bike$weekday, levels=0:6, labels = c('SUN', 'MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT'))
bike$holiday = factor(bike$holiday, levels = c(0,1), labels = c('NO HOLIDAY', 'HOLIDAY'))
bike$work = factor(bike$workingday, levels = c(0,1), labels = c('NO WORKING DAY', 'WORKING DAY'))
bike$workingday = NULL
bike$season = factor(bike$season, levels = 1:4, labels = c('WINTER', 'SPRING', 'SUMMER', 'FALL'))
bike$weather = factor(bike$weathersit, levels = 1:3, labels = c('GOOD', 'BAD', 'BAD'))
bike$weathersit = NULL
bike$mnth = factor(bike$mnth, levels = 1:12, labels = c('JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OKT', 'NOV', 'DEZ'))
bike$mnth = NULL
day_diff = function(date1, date2){
  as.numeric(difftime(as.Date(date1), as.Date(date2), units = 'days'))
}
bike$days_since_2011 = day_diff(bike$dteday, min(as.Date(bike$dteday)))
bike$yr[bike$yr == 0] = 2011
bike$yr[bike$yr == 1] = 2012
bike$yr = factor(bike$yr)
bike$registered = NULL
bike$casual = NULL
bike$instant = NULL
bike$dteday = NULL
bike$atemp = NULL
bike$temp = round(bike$temp * (39 - (-8)) + (-8), 1)
#windspeed: Normalized wind speed. The values are divided to 67 (max)
bike$wind = 67 * bike$windspeed
bike$windspeed = NULL
#hum: Normalized humidity. The values are divided to 100 (max)
bike$hum = round(100 * bike$hum, 3)
bike$days_since_2011 = NULL

#==============================================================================
# Measure predictive dependence
#==============================================================================
# Measure: Predictive dependence
feature_names = setdiff(colnames(bike), "cnt")
ind1 = sample(1:nrow(bike), size = nrow(bike) / 2)
ind2 = setdiff(1:nrow(bike), ind1)
r2s = lapply(feature_names, function(fname) {
  print(fname)
  form = as.formula(sprintf("%s ~ .", fname))
  r2 = c()
  for(ind in list(ind1, ind2)){
    mod = randomForest(form, data = bike[ind, feature_names])
    ind_other = setdiff(1:nrow(bike), ind)
    pred = predict(mod, newdata = bike[ind_other, feature_names])
    ytrue = bike[ind_other, fname]
    if (class(bike[[fname]]) == "numeric") {
      SSE = sum((pred - ytrue)^2)
      SST = sum((ytrue - mean(ytrue))^2)
    } else {
      tab = table(ytrue)
      mode_class = names(which.max(tab))
      mmce = function(y1, y2) sum(diag(table(y1,y2))) / length(y1)
      SST = 1 - mmce(ytrue, rep(mode_class, times = length(ytrue)))
      SSE = 1 - mmce(ytrue, pred)
    }
    r2 = c(r2, 1 - SSE/SST)
  }
  r2
})
r2sm = lapply(r2s, mean)
bike_pred_cor = data.frame(feature = feature_names, pred_cor = unlist(r2sm))
# Order for later plotting
bike_pred_cor$feature = factor(bike_pred_cor$feature, levels = bike_pred_cor$feature[order(bike_pred_cor$pred_cor)])
# Taking abs() of cor since for yr it was an epsilon below zero
bike_pred_cor$pred_cor = sprintf("%.0f%s", 100 * abs(bike_pred_cor$pred_cor), "%")
cap = "Percentage of loss explained by predicting a feature from the remaining features with a random forest."
lab = "tab:predcor"
xtab = xtable(t(bike_pred_cor), label = lab, caption = cap)
print(xtab, include.colnames = FALSE, include.rownames = FALSE, file = sprintf("%s/bike-pred-cor.tex", data_dir))

#==============================================================================
# Train model
#==============================================================================
# Split for training+tuning and for evaluation+iml
train_ids = sample(1:nrow(bike), size = (1 - test_size) * nrow(bike))
test_ids = setdiff(seq_len(nrow(bike)), train_ids)
task = TaskRegr$new(id = "bike", backend = bike, target = "cnt")
learner = lrn("regr.ranger")
mod = learner$train(task, row_ids = train_ids)
prediction = learner$predict(task, row_ids = test_ids)
perf = prediction$score(msr("regr.mae"))

#ctrl2 = rpart::rpart.control(maxdepth = 2)
feature_names = setdiff(colnames(bike), "cnt")
ctrl1 = partykit::ctree_control(maxdepth = 1, minbucket = 30)
conds1 = fit_conditionals(bike[train_ids, feature_names], ctrl = ctrl1, type = "trtf")
ctrl2 = partykit::ctree_control(maxdepth = 2, minbucket = 30)
conds2 = fit_conditionals(bike[train_ids, feature_names], ctrl = ctrl2, type = "trtf")
ctrl3 = partykit::ctree_control(maxdepth = 3, minbucket = 30)
conds3 = fit_conditionals(bike[train_ids, feature_names], ctrl = ctrl3, type = "trtf")
ctrl4 = partykit::ctree_control(maxdepth = 4, minbucket = 30)
conds4 = fit_conditionals(bike[train_ids, feature_names], ctrl = ctrl4, type = "trtf")
ctrl5 = partykit::ctree_control(maxdepth = 5, minbucket = 30)
cond5 = fit_conditionals(bike[train_ids, feature_names], ctrl = ctrl5, type = "trtf")


pred = Predictor$new(mod, data = bike[test_ids,], y = "cnt")

write(round(perf, 2), file = sprintf("%s/bike-mae.txt", data_dir))
write(round(cor(bike$temp, bike$hum), 2), file = sprintf("%s/bike-cor-hum-temp.txt", data_dir))



#==============================================================================
# subgroup PDP example
#==============================================================================

clean_path = function(pathstr){
  pathstr = gsub("[, ]?[\"]NA[\"][ ,]?", "", pathstr)
  pathstr = gsub("%in%", "in", pathstr, fixed = TRUE)
  pathstr = gsub("c(", "{", pathstr, fixed = TRUE)
  pathstr = gsub(")", "}", pathstr, fixed = TRUE)
  pathstr = gsub("\"", "", pathstr)
  pathstr = gsub("season in {WINTER, SPRING, FALL} &\n season in {WINTER,}", "season in {winter}", pathstr, fixed = TRUE)
  pathstr = gsub("season in {WINTER, SPRING, FALL} &\n season in {SPRING, FALL}", "season in {spring,fall}", pathstr, fixed = TRUE)
  pathstr = gsub("WORKING DAY", "WORK", pathstr, fixed = TRUE)
  pathstr = gsub("MON, TUE, WED, THU, FRI", "MO, ..., FRI", pathstr, fixed = TRUE)
  pathstr = gsub("workingday", "work", pathstr)
  pathstr = gsub("season in {SUMMER}", "season in {summer}", pathstr, fixed = TRUE)
  pathstr = gsub("season in {WINTER, SPRING, FALL} &\n season in {SPRING,}", "season in {spring}", pathstr, fixed = TRUE)
  pathstr = gsub("season in {WINTER, FALL} &\n season in {WINTER,}", "season in {winter}", pathstr, fixed = TRUE)
  pathstr = gsub("season in {WINTER, FALL} &\n season in { FALL}", "season in {fall}", pathstr, fixed = TRUE)
  pathstr = gsub("season in {SPRING, SUMMER} &\n season in { SUMMER,}", "season in {summer}", pathstr, fixed = TRUE)
  pathstr = gsub("season in {SPRING, SUMMER} &\n season in { SPRING,}", "season in {spring}", pathstr, fixed = TRUE)
  pathstr = gsub("season in {WINTER, SPRING, FALL} &\n season in {WINTER, FALL}", "season in {winter,fall}", pathstr, fixed = TRUE)
  pathstr
}


fname = "season"
cs_dat  = grouped_pdp_factor(pred, fname, conds2[[fname]])$df
cs_dat$node = clean_path(cs_dat$node)
p1 = ggplot(data.frame(cs_dat)) +
  geom_boxplot(aes_string(x = fname, y = ".value")) +
  facet_wrap("node", scales = "free_x", nrow = 1) +
  scale_y_continuous("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 = FeatureEffect$new(pred, fname, method = "ice")$plot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous("Predicted rentals")
cat_plot = p2 + p1 + plot_layout(widths = c(1,3))
pdf(file = sprintf("%s/bike-cat-plot.pdf", fig_dir), height = 4, width = 8)
print(cat_plot)
dev.off()



fname = "temp"
th = theme_bw() +
    theme(strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = "white",colour = "white"),
          text = element_text(size=18),
          legend.position = "top")


lmts = c(2000, 6500) 
yscale1 = scale_y_continuous("Predicted bike rentals", limits = lmts)
yscale2 = scale_y_continuous("", limits = lmts)
cs_dat  = grouped_pdp(pred, fname, conds2[[fname]])
cs_dat$pdp$node = clean_path(cs_dat$pdp$node)
cs_dat$outliers$node = clean_path(cs_dat$outliers$node)
p1 = ggplot(data = cs_dat$pdp,
           aes_string(x = fname, y = ".value",
                      group = "node", color = "node")) +
      geom_line(data = cs_dat$pdp[cs_dat$pdp$in_hinges,]) +
      geom_line(data = cs_dat$pdp[cs_dat$pdp$in_box,], size = 2) +
      geom_point(data = cs_dat$outliers) +
      yscale2 +
      scale_color_viridis("", discrete = TRUE) +
      th +
      scale_x_continuous("Temperature") +
      guides(color = guide_legend(nrow = 2))
tpdp = FeatureEffect$new(pred, fname, method = "pdp")$results
ale_results = FeatureEffect$new(pred, "temp", grid.size = 25)$results
mp = mean(pred$predict(bike)[[1]])
ale_results$.value = ale_results$.value + mp

effect = rbind(data.table(tpdp), data.table(ale_results))

effect_plot = ggplot(effect) +
  geom_line(aes_string(x = fname, y = ".value", group = ".type", lty = ".type")) +
  yscale2 +
  scale_linetype("") + th +
  scale_x_continuous("Temperature")
temp_plot = effect_plot + p1 +  plot_layout(widths = c(1, 2))
pdf(file = sprintf("%s/bike-temp-plot.pdf", fig_dir), height = 5, width = 12)
print(temp_plot)
dev.off()

#==============================================================================
# Feature Importance by depth
#==============================================================================
fimp = FeatureImp$new(pred, loss = "mae", n.repetitions = 50, compare = "difference")$results
fimp$depth = 0
fimp1 = grouped_pfi(pred, Metrics::mae, conds1, repetitions = 50)
fimp1$depth = 1
fimp2 = grouped_pfi(pred, Metrics::mae, conds2, repetitions = 50)
fimp2$depth = 2
fimp3 = grouped_pfi(pred, Metrics::mae, conds3, repetitions = 50)
fimp3$depth = 3
fimp4 = grouped_pfi(pred, Metrics::mae, conds4, repetitions = 50)
fimp4$depth = 4
fimp5 = grouped_pfi(pred, Metrics::mae, conds4, repetitions = 50)
fimp5$depth = 5

fimps = rbind(data.table(fimp), fimp1, fimp2, fimp3, fimp4, fimp5, fill = TRUE)

p = ggplot(fimps) +
  geom_line(aes(x = depth, y = importance, group = feature, color = feature)) +
  scale_y_continuous("cond. PFI") +
  scale_x_continuous("Maxdepth parameter")
filename = sprintf("%s/bike-importance-depth.pdf", fig_dir)
ggsave(file = filename, plot = p, width = 6, height = 3)
#==============================================================================
# Feature Importance Plot
#==============================================================================


fimp$type = "PFI"
cfimp = fimp2
cfimp$type = "cs-PFI"
fimps = rbind(data.table(fimp), cfimp, fill = TRUE)
fimps$feature = factor(fimps$feature, levels = rev(fimp$feature))
cfimp_g = grouped_pfi_groupwise(pred, Metrics::mae, conds2, repetitions = 50)
cfimp_g$feature_ylab = clean_path(cfimp_g$group)

fimps$condition[fimps$condition == ""] = NA


top5 = c("temp", "yr", "season", "hum", "wind")
fimps = fimps[fimps$feature %in% top5,]

p1 = ggplot(fimps, aes(x = importance, y = feature, color = type)) +
	geom_point() +
	geom_segment(aes(y = feature, yend = feature, x = importance.05, xend = importance.95)) +
  geom_text(aes(y = feature, x = importance, label = condition), show.legend = FALSE, nudge_y = -0.3, size = 3.5, hjust = 0) +
  scale_x_continuous("PFI and cs-PFI", limits = c(NA, max(fimps$importance) * 1.4)) +
  scale_y_discrete("") +
  scale_color_viridis("", discrete = TRUE) +
  theme(legend.position = "top")

p2 = ggplot(cfimp_g[cfimp_g$feature %in% c("temp"),], aes(x = importance, y = feature_ylab)) +
  geom_point() +
  geom_segment(aes(yend = feature_ylab, x = importance.05, xend = importance.95)) +
  scale_x_continuous("temperature cs-PFIs") +
  scale_y_discrete("")

pdf(file = sprintf("%s/bike-importance.pdf", fig_dir), height = 4, width = 8)
print(p1 + p2)
dev.off()


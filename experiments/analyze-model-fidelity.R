# =============================================================================
# Analyze data from the model fidelity experiment
# =============================================================================

devtools::load_all()
dataffiles = list.files(sprintf("%s/model-fidelity", res_dir), full.names = TRUE)
ex = lapply(dataffiles, read.csv)
ex = rbindlist(ex)

rnm_dat = c("space_ga" = "space",
            "wind" = "wind",
            "satellite_image" = "satellite",
            "pollen" = "pollen",
            "wine_quality" = "wine",
            "quake" = "quake")

rnm_levels = c("pdp_loss" = "PDP",
               "ale_loss" = "ALE",
               "gcpdp_loss1" = "cs-PDP trtr1",
               "gcpdp_loss1c"  = "cs-PDP cart1",
               "gcpdp_loss2" = "cs-PDP trtr2",
               "gcpdp_loss2c"  = "cs-PDP cart2",
               "gcpdp_loss5" = "cs-PDP trtr5",
               "gcpdp_loss5c"  = "cs-PDP cart5",
               "gcpdp_loss10" = "cs-PDP trtr10",
               "gcpdp_loss10c" = "cs-PDP cart10")

models = c("regr.ranger" = "RF",
           "regr.lm" = "LM",
           "regr.kknn" = "KN")


ex = melt(data.table(ex), measure.vars = names(rnm_levels))
# We only look at few groups here, 2 and 4
rmv = c("gcpdp_loss5", "gcpdp_loss10", "gcpdp_loss5c", "gcpdp_loss10c")
ex = ex[!(ex$variable %in% rmv), ]
ex$variable = rnm_levels[ex$variable]
ex$variable = factor(ex$variable, levels = c("PDP", "ALE", "cs-PDP cart1","cs-PDP trtr1",
                                             "cs-PDP cart2",  "cs-PDP trtr2"))


ex$learner = models[ex$learner]
ex$dat_name = rnm_dat[ex$dat_name]
pdf(file = sprintf("%s/model-fidelity.pdf", fig_dir), width = 6, height = 6)
p = ggplot(ex) +
	geom_boxplot(aes(x = variable, y = value, group = variable)) +
	facet_grid(dat_name ~ learner, scales = "free") +
        scale_x_discrete("") +
        scale_y_continuous("MSE") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(p)
dev.off()



# Create table with MSE. Only for random forest for space reasons.
ex2 = ex
ex2 = filter(ex, learner == "RF")
exs  = ex2 %>% group_by(dat_name, learner, variable) %>%
  summarize(median(value))

exs_wide = dcast(data.table(exs), variable ~ dat_name)
method = exs_wide[,1][[1]]
exs_wide[,1] = NULL
exs_wide = round(exs_wide, 4)
exs_wide = data.frame(exs_wide)
rownames(exs_wide) = method
# Sort by method
exs_wide = exs_wide[rnm_levels,]
exs_wide = na.omit(exs_wide)
cap = "Median model fidelity averaged over features in a random forest for various data sets. The cPDPs always had a lower loss (i.e. higher model fidelity) than PDP and ALE. The loss monotonically decreases with increasing maximum tree depth for subgroup construction."
lab = "tab:model-fidelity"
tab = xtable(exs_wide, format = "latex", caption = cap, label = lab)
print(tab, file = sprintf("%s/model-fidelity-tab.tex", fig_dir))


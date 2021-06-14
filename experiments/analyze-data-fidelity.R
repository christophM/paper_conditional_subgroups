# =============================================================================
# Analyze data from the data fidelity (MMD) experiment
# =============================================================================
devtools::load_all()

FIG_HEIGHT = 8
TEXT_SIZE = 9


selection_global = c("perm", "cvirf", "cart30", "trtr30", "ale", "ko", "imp", "none")
selection_tree = c("perm",
									 "trtr1", "cart1", "trtr2", "cart2", "trtr3", "cart3",
                   "trtr4", "cart4", "trtr5", "cart5", "trtr30", "cart30",
									 "none")

dataffiles = list.files(sprintf("%s/data-fidelity", res_dir), full.names = TRUE)
results = lapply(dataffiles, read.csv)
res = rbindlist(results, fill = TRUE)
colnames(res) = gsub("trtf", "trtr", colnames(res))
colnames(res) = gsub("strobl", "cvirf", colnames(res))
variables = unique(c(selection_global, selection_tree))


res$imp = res$imputation
res$imputation = NULL

res.m = melt(res, measure.var = variables)
res.m$variable = as.character(res.m$variable)

ranks = res.m %>%
  filter(variable %in% selection_global) %>%
	#filter(grid_type == "quantile") %>%
  group_by(data_name, feature, repetition) %>%
  mutate(rankz = rank(value)) %>%
  ungroup() %>%
	# Average over repetitions
	group_by(variable, data_name, feature) %>%
	mutate(rankz = mean(rankz)) %>%
	ungroup() %>%
	# Summarize over data sets and features
  group_by(variable) %>%
  summarize(mrank = mean(rankz), srank = sd(rankz))
ranks$mrank = round(ranks$mrank, 2)
ord = order(ranks$mrank)
tab = t(ranks[c("mrank", "srank")])
colnames(tab) = ranks$variable
rownames(tab) = c("Mean ranks", "SD")
tab = tab[,ord, drop = FALSE]
colnames(tab)[colnames(tab) == "cart30"] = "cs (cart)"
colnames(tab)[colnames(tab) == "trtr30"] = "cs (trtr)"

xtab = xtable(tab, label = "tab:ranks", caption = "Mean ranks and their standard deviation based on data fidelity of various perturbation methods over data sets, features and repetitions. \\textbf{Legend:} none: No intervention, which serves as upper benchmark. cart30: cs-permutation with CART with maximal depth of 30. trtr30: cs-permutation with transformation trees with maximal depth of 30. imp: Imputation approach. ko: Model-X knockoffs \\cite{candes2018panning} . ale: ALE perturbation \\cite{apley2016visualizing}. cvirf: Conditional variable importance for random forests \\cite{strobl2008conditional}. perm: Unconditional permutation.")
print(xtab, file = sprintf("%s/mmd-ranks.tex", fig_dir))


plot_data_fidelity = function(dat, ncol = 5) {
  ggplot(dat) +
    geom_boxplot(aes(x = variable, y = -log(value))) +
    xlab("") +
    ylab("Data Fidelity (-log(MMD)") +
    facet_wrap("data_name", scales = "free_y", ncol = ncol) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0),
          text = element_text(size = TEXT_SIZE))
}

compare = res.m %>%
  group_by(data_name, feature, variable) %>%
  summarize(value = mean(value)) %>%
  group_by(data_name, feature) %>%
  mutate(value = (-log(value) - min(-log(value))) / (max(-log(value)) - min(-log(value))))

compare = filter(compare, variable %in% selection_global)
ord = compare %>%
        group_by(variable) %>%
        summarize(value=mean(value)) %>%
        arrange(-value) %>% pull(variable)
compare$variable = factor(compare$variable, levels = intersect(ord, selection_global))

pdf(file = sprintf("%s/data-mmd-all.pdf", fig_dir), height = 5)
ggplot(compare) +
  geom_boxplot(aes(x = variable, y = value)) +
  xlab("") +
  ylab("scaled Data Fidelity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0),
        text = element_text(size = TEXT_SIZE))
dev.off()


# split dataset set in half two split plots
data_names = unique(res.m$data_name)
data1 = data_names[1:round(length(data_names)/2)]
data2 = setdiff(data_names, data1)

## Dataset-wise comparison of approaches
# Use only deepest trees
res.m_all = filter(res.m, variable %in% selection_global)
res.m_all$variable = factor(res.m_all$variable, levels = selection_global)
res.m_all  = res.m_all %>%
  group_by(data_name, feature, variable) %>%
  summarize(value = median(value))



pdf(file = sprintf("%s/data-mmd-all1.pdf", fig_dir), height = FIG_HEIGHT)
plot_data_fidelity(res.m_all[res.m_all$data_name %in% data1,], ncol = 6)
dev.off()

pdf(file = sprintf("%s/data-mmd-all2.pdf", fig_dir), height = FIG_HEIGHT)
plot_data_fidelity(res.m_all[res.m_all$data_name %in% data2,], ncol = 6)
dev.off()




plot_data_fidelity2 = function(dat, ncol = 5) {
  dat1 = dat[!(dat$variable %in% c("none", "perm")),]
  median_perm = dat %>%
    filter(variable == "perm") %>%
    group_by(data_name) %>%
    summarize(med = median(-log(value)))
  median_none = dat %>%
    filter(variable == "none") %>%
    group_by(data_name) %>%
    summarize(med = median(-log(value)))

  p = plot_data_fidelity(dat1, ncol = ncol)
  p + geom_hline(aes(yintercept = med), data = median_perm, lty = 2) +
    geom_hline(aes(yintercept = med), data = median_none, lty = 2)
}


## Dataset-wise comparison of CART vs. TRTF and different depths
res.m_tree = filter(res.m, variable %in% selection_tree)
res.m_tree$variable = factor(res.m_tree$variable, levels = selection_tree)
res.m_tree= res.m_tree %>%
  group_by(data_name, feature, variable) %>%
  summarize(value = median(value))

pdf(file = sprintf("%s/data-mmd-tree1.pdf", fig_dir), height = FIG_HEIGHT)
plot_data_fidelity2(res.m_tree[res.m_tree$data_name %in% data1,])
dev.off()

pdf(file = sprintf("%s/data-mmd-tree2.pdf", fig_dir), height = FIG_HEIGHT)
plot_data_fidelity2(res.m_tree[res.m_tree$data_name %in% data2,])
dev.off()



plot_lmer_cis = function(lsmeans, var_order = NULL, negative = FALSE){
  lsmeans$variable = gsub("variable", "", rownames(lsmeans))
  if (!is.null(var_order)){
    lsmeans$variable = factor(lsmeans$variable, levels = var_order)
  }
  if (negative) {
    lsmeans$Estimate = -lsmeans$Estimate
    lsmeans$lower = - lsmeans$lower
    lsmeans$upper = - lsmeans$upper
  }
  ggplot(lsmeans) +
    geom_segment(aes(y = variable, yend = variable, x = lower, xend = upper)) +
    geom_point(aes(x = Estimate, y = variable)) +
    scale_y_discrete("Perturbation method") +
    scale_x_continuous("Data fidelity estimates") +
    theme(text = element_text(size = 13))
}


# Mixed effect model
res2 = filter(res.m, variable %in% selection_global)
#res2 = filter(res2, !(variable == "ale"))
res2$variable = factor(res2$variable)
res2$variable = relevel(res2$variable, "perm")
mod = lmerTest::lmer(-log(value) ~ variable + (1 | data_name / feature), data = res2)
pairwise = ls_means(mod, pairwise = TRUE)
rownames(pairwise) = gsub("variable", "", rownames(pairwise))

rows = c("perm - cart30", "perm - trtr30", "perm - ale", "perm - cvirf", "perm - imp", "perm - ko", "perm - none")
lsmeans = pairwise[rows,]
rownames(lsmeans) = gsub("perm - ", "", rownames(lsmeans))
# Order by descending data fidelity
var_order = summary(mod)$coefficients[,"Estimate"]
var_order = c(var_order, "perm"=0)
var_order = names(var_order)[order(var_order)]
var_order = setdiff(var_order, "(Intercept)")
var_order = gsub("variable", "", var_order)
var_order[var_order == "cart30"] = "cs (cart)"
var_order[var_order == "trtr30"] = "cs (trtr)"

rownames(lsmeans)[rownames(lsmeans) == "cart30"] = "cs (cart)"
rownames(lsmeans)[rownames(lsmeans) == "trtr30"] = "cs (trtr)"



# Compare QQ norm
qqnorm(residuals(mod))
qqline(residuals(mod))

fig_a =  plot_lmer_cis(lsmeans, var_order = var_order, negative = TRUE)
fig_a = fig_a + ggtitle("A)")

res3 = filter(res.m, variable %in% selection_tree)
res3$variable = factor(res3$variable)
res3$variable = relevel(res3$variable, "perm")

mod = lmerTest::lmer(-log(value) ~  variable  + (1 | data_name / feature), data = res3)
## Pairwise comparisons for different depths
lsmeans = ls_means(mod, pairwise = TRUE)
rownames(lsmeans) = gsub("variable", "", rownames(lsmeans))
rows = c("perm - trtr1", "perm - cart1", "perm - cart2", "perm - trtr2", "perm - cart3", "perm - trtr3", "perm - cart4", "perm - trtr4", "perm - cart5", "perm - trtr5", "perm - cart30", "perm - trtr30", "perm - none")
lsmeans = lsmeans[rows,]
rownames(lsmeans) = gsub("perm - ", "", rownames(lsmeans))

# Order by descending data fidelity
var_order = summary(mod)$coefficients[,"Estimate"]
var_order = c(var_order, "perm"=0)
var_order = names(var_order)[order(var_order)]
var_order = setdiff(var_order, "(Intercept)")
var_order = gsub("variable", "", var_order)

# Compare QQ norm
qqnorm(residuals(mod))
qqline(residuals(mod))

fig_b = plot_lmer_cis(lsmeans, var_order = var_order, negative = TRUE)
fig_b = fig_b + ggtitle("B)") + scale_x_continuous(limits = c(0, NA))

pdf(file = sprintf("%s/data-fidelity-cis.pdf", fig_dir), height = 3)
print(fig_a + fig_b)
dev.off()






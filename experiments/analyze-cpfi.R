# =============================================================================
# Analyze data from the cPFI groundtruth experiment
# =============================================================================
devtools::load_all()



mlegend = "Legend: impute rf: Imputation with a random forest, ko: Model-X knockoffs, mPFI: (marginal) PFI, tree cart: cs-permutation based on CART, tree trtr: cs-permutation based on transformation trees"
mlegendplus = sprintf("%s, CVIRF: conditional variable importance for random forests", mlegend)

# Experiments comparing
ex = readRDS(sprintf("%s/true-experiment.RDS", res_dir))
res = lapply(ex$result, function(x){
  ifelse(length(x) == 0, NA, x[[1]])
})
ex$result = unlist(res)

# Check completeness of results
# ex %>%
#  group_by(problem, algorithm, p, n, g) %>%
# summarize(n()) %>% data.frame

ex2 = readRDS(sprintf("%s/true-experiment-depth.RDS", res_dir))
ex2$result = ex2$pfi
ex = rbind(ex, ex2, fill = TRUE)
ex$setting = sprintf("n=%i, p=%i",ex$n, ex$p)
ex$g = gsub("_g", "", ex$g)
ex$g = gsub("nlinear", "non-linear", ex$g)
ex$g = gsub("high_dim", "multi. lin.", ex$g)
ex$type = as.character(ex$type)
ex$type[ex$type == "trtf"] = "trtr"

ex = filter(ex, !(algorithm == "imp_lm"))

# =============================================================================
# Merge true PFI (true)
# =============================================================================
true_pfi = dplyr::filter(ex, algorithm == "true", problem == "e1") %>%
  dplyr::group_by(problem, f, g, setting, n, p) %>%
  dplyr::summarize(tpfi = mean(result), .groups = "drop") 
ex = filter(ex, !(algorithm == "true"))
ex = merge(ex, true_pfi, by = c("problem", "g", "f", "setting", "n", "p"), all.x = TRUE)

# =============================================================================
# Merge true of the random forest model (true_rf)
# =============================================================================

# Uncomment me if experiment-cpfi.R is rerun
#true_pfi_rf = dplyr::filter(ex, algorithm == "true", problem == "e2") %>%
#  dplyr::group_by(problem, f, g, setting, n, p) %>%
#  dplyr::summarize(tpfi_rf = mean(result), .groups = "drop") 
#ex = filter(ex, !(algorithm == "true_rf"))
#ex = merge(ex, true_pfi_rf, by = c("problem", "g", "f", "setting", "n", "p"), all.x = TRUE)

# Delete me if experiment-cpfi.R is rerun
# <<<
## Code to produce csv for the extra "manual" run of experiment-cpfi where the bug in true_rf was fixed
filename = sprintf("%s/results/true-rf.csv", here())
#true_rf = filter(ex, problem == "e2", algorithm == "true_rf")
#true_rf %>% group_by(problem, f, p, n, g) %>%
#  summarize(tpfi_rf = mean(result), .groups = "drop") %>%
#  arrange(g, p, n) %>%
#  data.frame() %>%
#  write.csv(file = filename)
true_pfi_rf = read.csv(filename)
true_pfi_rf$setting = sprintf("n=%i, p=%i", true_pfi_rf$n, true_pfi_rf$p)
ex = filter(ex, !(algorithm == "true_rf"))
ex = merge(ex, true_pfi_rf, by = c("problem", "g", "f", "n", "p", "setting"), all.x = TRUE)
# delete until here
# >>>

# =============================================================================
# General settings
# =============================================================================
th = theme(text = element_text(size=18), legend.position = "top")

# Rename the algorithms
ex$algorithm = as.character(ex$algorithm)
rnm = c("tree" = "tree", "imp" = "impute rf", "imp_lm" = "impute lm",
        "ko" = "ko", "marg" = "mPFI", "cvirf" = "cvirf")
ex$algorithm = rnm[ex$algorithm]

# =============================================================================
# Tree Depth Experiments
# =============================================================================
ex1 = filter(ex, algorithm == "tree", problem == "e1", !is.na(nterminal))
ex_summary = ex1 %>% 
  filter(minbucket == 30, problem == "e1", nterminal > 0) %>%
  group_by(type, setting, p, g,n, nterminal, minbucket) %>%
  filter(n() > 30) %>%
  summarize(n = n(),
            pfi_sd = sd(pfi)/n(), size = n(),
            lower = quantile(pfi, prob = 0.05),
            upper = quantile(pfi, prob = 0.95),
            pfi = median(pfi),
  )

p = ggplot(aes(x = nterminal, y = pfi,  group = type), data = ex_summary) +
  geom_line(aes(color = type, lty = type), size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = type), alpha = 0.3) +
  geom_hline(aes(yintercept = tpfi), data = true_pfi[true_pfi$problem == "e1",]) +
  scale_y_continuous("PFI") +
  scale_x_continuous("Number of subgroups") +
  scale_color_discrete("Tree algorithm") +
  scale_linetype_discrete("Tree algorithm") +
  scale_fill_discrete("Tree algorithm") +
  th +
  guides(color = guide_legend(nrow = 2)) +
  facet_grid(g ~ setting, scales = "free")

print(p)
ggsave(file = sprintf("%s/true-importance-trees-mb30.pdf", fig_dir), plot = p, width = 12, height = 12)

# =============================================================================
# Compare for true f 
# =============================================================================

# ggsave
comp_all = filter(ex, (algorithm != "tree") | ((maxdepth == 30) & (minbucket == 30)))
comp_all$algorithm[comp_all$algorithm == "tree"] = sprintf("cs-PFI (%s)",comp_all$type[comp_all$algorithm == "tree"])

comp_all1 = filter(comp_all, problem == "e1", algorithm != "cvirf")
p2_true = ggplot(aes(x = algorithm, y = result, fill = algorithm), data = comp_all1) +
  geom_violin() +
  geom_hline(aes(yintercept = tpfi), data = true_pfi) +
  scale_fill_viridis("", discrete = TRUE) +
  scale_x_discrete("") +
  scale_y_continuous("PFI") +
  guides(fill = guide_legend(nrow = 1)) +
  facet_grid(g ~ setting, scales = "free_y") +
  th +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
  legend.position = "none")
# ggsave

ggsave(file = sprintf("%s/true-importance-all-ex1.pdf", fig_dir), plot = p2_true, width = 12, height = 12)
print(p2_true)


# Table format
tab = comp_all1 %>% group_by(g, setting, algorithm) %>%
  summarize(mse = mean((result- tpfi)^2)) %>%
  tidyr::pivot_wider(names_from = algorithm, values_from = mse)
cap = "MSE comparing estimated and true conditional PFI (scenario I)."
cap = sprintf("%s %s.", cap, mlegend)
lab = "mses-ex1"
filename = sprintf("%s/mses-ex1.tex", fig_dir)
pack_index = c(4, 4, 4, 4)
names(pack_index) = unique(tab$g)
tab$g = NULL
kbl(tab, booktabs = TRUE, format = "latex", digits = 2,
      caption = cap, label = lab, escape = FALSE, linesep = "") %>%
  kable_classic() %>%
  pack_rows(index = pack_index) %>%
  write(file = filename)

# =============================================================================
# Compare for intermediate conditional forest
# =============================================================================
comp_all2 = filter(comp_all, problem == "e2", is.na(permimp_threshold) | (permimp_threshold == 0.95))
p2_true_rf = ggplot(aes(x = algorithm, y = result, fill = algorithm), data = comp_all2) +
  geom_violin() +
  geom_hline(aes(yintercept = tpfi_rf), data = true_pfi_rf) +
  scale_fill_viridis("", discrete = TRUE) +
  scale_x_discrete("") +
  scale_y_continuous("PFI") +
  guides(fill = guide_legend(nrow = 1)) +
  facet_grid(g ~ setting, scales = "free_y") +
  th +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
  legend.position = "none")

# ggsave
plot(p2_true_rf)
ggsave(file = sprintf("%s/true-importance-all-ex2.pdf", fig_dir), plot = p2_true_rf, width = 12, height = 10)

# Table format
tab = comp_all2 %>% group_by(g, setting, algorithm) %>%
  summarize(mse = mean((result - tpfi_rf)^2)) %>%
  tidyr::pivot_wider(names_from = algorithm, values_from = mse)
cap = "MSE comparing estimated and true conditional PFI (for random forest, scenario II)."
cap = sprintf("%s %s.", cap, mlegendplus)
lab = "mses-ex2"
filename = sprintf("%s/mses-ex2.tex", fig_dir)
pack_index = c(4, 4, 4, 4)
names(pack_index) = unique(tab$g)
tab$g = NULL
kbl(tab, booktabs = TRUE, format = "latex", digits = 2,
      caption = cap, label = lab, escape = FALSE, linesep = "") %>%
  kable_classic() %>%
  pack_rows(index = pack_index) %>%
  write(file = filename)


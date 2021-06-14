# =============================================================================
# Illustrate boxplot PDP variant 
# =============================================================================

suppressPackageStartupMessages({devtools::load_all()})
set.seed(42)
# caption: The evaluation of each point at the PDP is the same, only the emphasis changes.
n = 100
x = rgamma(n, 1)
fy = function(x) log(x) + 4.5
d = data.frame(x = x, y = fy(x))


p1 = ggplot(d) + geom_line(aes(x = x, y = y))

bpstats = data.frame(stats = boxplot.stats(x)$stats)
bpstats$y = fy(bpstats$stats) 
outliers = x[x < bpstats$stats[1] |
             x > bpstats$stats[5]]
outliers = data.frame(x = outliers, y = fy(outliers))
in_hinges =  x >= bpstats$stats[1] &
             x <= bpstats$stats[5]
in_box = x >= bpstats$stats[2] &
         x <= bpstats$stats[4]


p2 = ggplot(mapping = aes(x = x, y = y)) + 
  geom_line(data = d[in_hinges,]) + 
  geom_line(data = d[in_box,], size = 2) + 
  geom_point(data = outliers) +
  scale_y_continuous("")

p3 = ggplot(d) + geom_boxplot(aes(x = x)) + theme_void()

layout = "
AAAAAABBBBBB
AAAAAABBBBBB
AAAAAABBBBBB
AAAAAABBBBBB
######CCCCCC
"
p = p1 + p2 + p3 + plot_layout(design = layout)
pdf(file = sprintf("%s/pdp-boxplot.pdf", fig_dir), height = 2)
print(p)
dev.off()



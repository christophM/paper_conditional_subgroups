# =============================================================================
# Global settings for R scripts
# =============================================================================
library("viridis")
library("ggplot2")

data_dir = here("./data/")
fig_dir = here("./paper/figures")
res_dir = here("./results")

# define color
gg_color_hue = function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
hue = gg_color_hue(2)

## ggplot black and white
theme_set(theme_bw() +
    theme(strip.text = element_text(face = "bold"),
		   strip.background = element_rect(fill = "white",colour = "white"),
		   text = element_text(size=14)))


# Overall controls for ctree
ctrl = function(md) partykit::ctree_control(maxdepth = md, minbucket = 30)
